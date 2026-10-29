use bio::bio_types::strand::Strand;

use crate::{
    constants::{
        CODON_LENGTH, MINIMUM_EDGE_GENE_LENGTH, MINIMUM_GENE_LENGTH, READING_FRAMES, STT_NOD,
    },
    node::validation::{
        check_start_codon, is_edge_gene, is_reverse_edge_gene, is_valid_gene, is_valid_reverse_gene,
    },
    sequence::{encoded::EncodedSequence, is_stop},
    types::{CodonType, Mask, Node, NodePosition, OrphosError, Training},
};

/// Type alias for arrays indexed by reading frame
type ReadingFrameArray<T> = [T; READING_FRAMES];

/// Context for processing DNA strands containing all tracking state
struct StrandProcessingContext {
    last_stop_positions: ReadingFrameArray<usize>,
    has_start_codon: ReadingFrameArray<bool>,
    minimum_distances: ReadingFrameArray<usize>,
    sequence_length: usize,
    closed: bool,
    circular: bool,
}

impl StrandProcessingContext {
    /// Create a new processing context with initialized arrays
    fn new(
        sequence_length: usize,
        sequence_length_mod: usize,
        closed: bool,
        circular: bool,
    ) -> Self {
        let mut context = Self {
            last_stop_positions: ReadingFrameArray::default(),
            has_start_codon: ReadingFrameArray::default(),
            minimum_distances: ReadingFrameArray::default(),
            sequence_length,
            closed,
            circular,
        };

        context.initialize_arrays(sequence_length_mod);
        context
    }

    /// Initialize the tracking arrays for strand processing
    fn initialize_arrays(&mut self, sequence_length_mod: usize) {
        // Initialize arrays with different indexing patterns to match original C code
        for i in 0..READING_FRAMES {
            // Different indexing patterns to match C code exactly
            self.last_stop_positions[(i + sequence_length_mod) % READING_FRAMES] =
                self.sequence_length + i; // Uses (i+slmod)%3 for last
            self.has_start_codon[i % READING_FRAMES] = false; // Uses i%3 for saw_start  
            self.minimum_distances[i % READING_FRAMES] = MINIMUM_EDGE_GENE_LENGTH; // Uses i%3 for min_dist

            if !self.closed && self.sequence_length > 0 {
                while self.last_stop_positions[(i + sequence_length_mod) % READING_FRAMES] + 2
                    > self.sequence_length - 1
                {
                    self.last_stop_positions[(i + sequence_length_mod) % READING_FRAMES] = self
                        .last_stop_positions[(i + sequence_length_mod) % READING_FRAMES]
                        .saturating_sub(3);
                }
            }
        }
    }

    fn seed_circular_stops(&mut self, encoded_sequence: &[u8], training: &Training) {
        if !self.circular || self.sequence_length < CODON_LENGTH {
            return;
        }

        let mut first_stop_per_frame: [Option<usize>; READING_FRAMES] = [None; READING_FRAMES];
        let scanning_start_position = self.sequence_length - CODON_LENGTH;

        for position_index in 0..=scanning_start_position {
            if is_stop(encoded_sequence, position_index, training) {
                let frame = position_index % READING_FRAMES;
                if first_stop_per_frame[frame].is_none() {
                    first_stop_per_frame[frame] = Some(position_index);
                }
            }
        }

        for (frame, first_stop) in first_stop_per_frame.iter().enumerate() {
            if let Some(stop_position) = first_stop {
                self.last_stop_positions[frame] = *stop_position;
                self.minimum_distances[frame] = MINIMUM_GENE_LENGTH;
            }
        }
    }
}

/// Classification of every 6-bit codon code for one translation table.
///
/// The tables are derived from [`is_stop`] and [`check_start_codon`] on each of
/// the 64 possible codons, so the strand scans classify codons exactly like
/// the per-position checks while reading every base only once.
struct CodonTables {
    is_stop: [bool; 64],
    start_type: [Option<CodonType>; 64],
}

impl CodonTables {
    fn new(training: &Training) -> Self {
        let mut tables = Self {
            is_stop: [false; 64],
            start_type: [None; 64],
        };
        for code in 0..64u8 {
            let (first, second, third) = (code >> 4, (code >> 2) & 0b11, code & 0b11);
            // Pack the three bases as the encoded sequence stores them.
            let packed = [first | (second << 2) | (third << 4), 0];
            tables.is_stop[usize::from(code)] = is_stop(&packed, 0, training);
            tables.start_type[usize::from(code)] = check_start_codon(&packed, 0, training);
        }
        tables
    }
}

/// Classification of the four codons that start inside one byte of the
/// encoded sequence, keyed by that byte and the low half of the following
/// one (six bases).
///
/// Entry layout: bits 0..4 flag stop codons starting at base 0..3 of the
/// byte, bits 4..8 flag start codons, and bits 8..16 hold the start codon
/// type index (two bits per base) for the flagged starts.
struct ByteCodonTable {
    entries: Vec<u16>,
}

impl ByteCodonTable {
    /// Six bases of two bits each.
    const KEY_BITS: usize = 12;
    const KEY_COUNT: usize = 1 << Self::KEY_BITS;
    const STOP_SHIFT: u16 = 0;
    const START_SHIFT: u16 = 4;
    const TYPE_SHIFT: u16 = 8;

    fn new(codon_tables: &CodonTables) -> Self {
        let entries = (0..Self::KEY_COUNT as u16)
            .map(|key| {
                let base = |offset: u16| ((key >> (2 * offset)) & 0b11) as u8;
                let mut entry = 0u16;
                for offset in 0..4u16 {
                    let codon = usize::from(
                        (base(offset) << 4) | (base(offset + 1) << 2) | base(offset + 2),
                    );
                    if codon_tables.is_stop[codon] {
                        entry |= 1 << (Self::STOP_SHIFT + offset);
                    }
                    if let Some(codon_type) = codon_tables.start_type[codon] {
                        entry |= 1 << (Self::START_SHIFT + offset);
                        entry |= (codon_type.to_index() as u16) << (Self::TYPE_SHIFT + 2 * offset);
                    }
                }
                entry
            })
            .collect();
        Self { entries }
    }

    /// Entry for the byte at `byte_index`; bases past the end of the buffer
    /// read as zero, which only affects codons past the scanned range.
    #[inline]
    fn entry(&self, encoded_sequence: &[u8], byte_index: usize) -> u16 {
        let low = usize::from(encoded_sequence[byte_index]);
        let high = usize::from(encoded_sequence.get(byte_index + 1).copied().unwrap_or(0));
        self.entries[(low | (high << 8)) & (Self::KEY_COUNT - 1)]
    }

    #[inline]
    fn is_stop(entry: u16, offset: usize) -> bool {
        entry & (1 << (Self::STOP_SHIFT as usize + offset)) != 0
    }

    #[inline]
    fn start_type(entry: u16, offset: usize) -> Option<CodonType> {
        if entry & (1 << (Self::START_SHIFT as usize + offset)) == 0 {
            return None;
        }
        Some(
            match (entry >> (Self::TYPE_SHIFT as usize + 2 * offset)) & 0b11 {
                0 => CodonType::Atg,
                1 => CodonType::Gtg,
                _ => CodonType::Ttg,
            },
        )
    }
}

/// Number of sequence positions encoded per byte.
const POSITIONS_PER_BYTE: usize = 4;

/// Add nodes for start and stop codons in both directions
///
/// This function identifies potential start and stop codons in the sequence
/// and creates corresponding nodes for gene prediction analysis.
///
/// # Arguments
/// * `encoded_sequence` - The complete encoded sequence data
/// * `nodes` - Vector to store the created nodes
/// * `closed` - Whether to treat sequence as circular/closed
/// * `circular` - Whether to treat sequence topology as circular for wraparound genes
/// * `training` - Training data for scoring parameters
///
/// # Returns
/// Number of nodes created, or error if processing fails
pub fn add_nodes(
    encoded_sequence: &EncodedSequence,
    nodes: &mut Vec<Node>,
    closed: bool,
    circular: bool,
    training: &Training,
) -> Result<usize, OrphosError> {
    let sequence_length = encoded_sequence.sequence_length;
    // Clear the nodes vector
    nodes.clear();
    // Enough for the start and stop codons of a typical genome (one node
    // every twenty bases or so), so the vector rarely has to grow and copy.
    nodes.reserve(STT_NOD.max(sequence_length / 16));

    let sequence_length_mod = sequence_length % READING_FRAMES;

    let codon_tables = ByteCodonTable::new(&CodonTables::new(training));

    // The two strands are independent; their nodes are created on two
    // threads, the forward strand's first as before. Sorting follows, so
    // the order within the vector does not matter beyond that.
    let (_, mut reverse_nodes) = rayon::join(
        || {
            let mut context = StrandProcessingContext::new(
                sequence_length,
                sequence_length_mod,
                closed,
                circular,
            );
            context.seed_circular_stops(&encoded_sequence.forward_sequence, training);
            process_forward_strand(
                &encoded_sequence.forward_sequence,
                nodes,
                &mut context,
                &encoded_sequence.masks,
                training,
                &codon_tables,
            );
            handle_remaining_starts(
                &encoded_sequence.forward_sequence,
                nodes,
                &context,
                Strand::Forward,
                training,
            );
        },
        || {
            let mut reverse_nodes = Vec::with_capacity(STT_NOD.max(sequence_length / 32));
            let mut context = StrandProcessingContext::new(
                sequence_length,
                sequence_length_mod,
                closed,
                circular,
            );
            context.seed_circular_stops(&encoded_sequence.reverse_complement_sequence, training);
            process_reverse_strand(
                &encoded_sequence.reverse_complement_sequence,
                &mut reverse_nodes,
                &mut context,
                &encoded_sequence.masks,
                training,
                &codon_tables,
            );
            handle_remaining_starts(
                &encoded_sequence.reverse_complement_sequence,
                &mut reverse_nodes,
                &context,
                Strand::Reverse,
                training,
            );
            reverse_nodes
        },
    );
    nodes.append(&mut reverse_nodes);

    Ok(nodes.len())
}

/// Process the forward strand to identify potential genes
///
/// Scans the forward strand from 3' to 5' end, looking for start and stop codons
/// and creating nodes for valid gene boundaries that meet length requirements.
fn process_forward_strand(
    encoded_sequence: &[u8],
    nodes: &mut Vec<Node>,
    context: &mut StrandProcessingContext,
    masks: &[Mask],
    training: &Training,
    codon_tables: &ByteCodonTable,
) {
    if context.sequence_length < CODON_LENGTH {
        return;
    }
    let scanning_start_position = context.sequence_length - CODON_LENGTH;

    // Scan one byte (four positions) at a time, from the end of the sequence.
    // Bytes without stop or start codons only matter for edge genes at the
    // first positions, so everything else is skipped without per-position
    // work.
    let mut byte_index = scanning_start_position / POSITIONS_PER_BYTE;
    loop {
        let entry = codon_tables.entry(encoded_sequence, byte_index);
        let first_position = byte_index * POSITIONS_PER_BYTE;
        let last_offset = (scanning_start_position - first_position).min(POSITIONS_PER_BYTE - 1);
        if entry != 0 || first_position <= EDGE_POSITION_LIMIT {
            for offset in (0..=last_offset).rev() {
                let position_index = first_position + offset;
                process_forward_position(
                    encoded_sequence,
                    nodes,
                    context,
                    masks,
                    training,
                    position_index,
                    ByteCodonTable::is_stop(entry, offset),
                    ByteCodonTable::start_type(entry, offset),
                );
            }
        }
        if byte_index == 0 {
            break;
        }
        byte_index -= 1;
    }
}

/// Positions at or below this may start an edge gene.
const EDGE_POSITION_LIMIT: usize = 2;

#[inline]
#[allow(clippy::too_many_arguments)]
fn process_forward_position(
    encoded_sequence: &[u8],
    nodes: &mut Vec<Node>,
    context: &mut StrandProcessingContext,
    masks: &[Mask],
    training: &Training,
    position_index: usize,
    is_stop_codon: bool,
    start_type: Option<CodonType>,
) {
    {
        let reading_frame_index = position_index % READING_FRAMES;

        if is_stop_codon {
            if context.has_start_codon[reading_frame_index] {
                let node = create_stop_node(
                    context.last_stop_positions[reading_frame_index],
                    position_index as isize,
                    Strand::Forward,
                    encoded_sequence,
                    training,
                );
                nodes.push(node);
            }

            context.minimum_distances[reading_frame_index] = MINIMUM_GENE_LENGTH;
            context.last_stop_positions[reading_frame_index] = position_index;
            context.has_start_codon[reading_frame_index] = false;
        } else if context.last_stop_positions[reading_frame_index] < context.sequence_length {
            // A start codon that fails the length or mask checks still falls
            // through to the edge-gene test, as in Prodigal.
            let valid_start = start_type.filter(|_| {
                is_valid_gene(
                    position_index,
                    context.last_stop_positions[reading_frame_index],
                    context.minimum_distances[reading_frame_index],
                    context.sequence_length,
                    context.circular,
                    masks,
                )
            });
            if let Some(codon_type) = valid_start {
                let node = create_start_node(
                    position_index,
                    codon_type,
                    context.last_stop_positions[reading_frame_index] as isize,
                    Strand::Forward,
                );
                context.has_start_codon[reading_frame_index] = true;
                nodes.push(node);
            } else if is_edge_gene(
                position_index,
                context.last_stop_positions[reading_frame_index],
                context.closed || context.circular,
                masks,
            ) {
                let node = create_edge_node(
                    position_index,
                    context.last_stop_positions[reading_frame_index] as isize,
                    Strand::Forward,
                );
                context.has_start_codon[reading_frame_index] = true;
                nodes.push(node);
            }
        }
    }
}

/// Process reverse strand to find genes
fn process_reverse_strand(
    reverse_complement_encoded_sequence: &[u8],
    nodes: &mut Vec<Node>,
    context: &mut StrandProcessingContext,
    masks: &[Mask],
    training: &Training,
    codon_tables: &ByteCodonTable,
) {
    if context.sequence_length < CODON_LENGTH {
        return;
    }
    let scanning_start_position = context.sequence_length - CODON_LENGTH;

    let mut byte_index = scanning_start_position / POSITIONS_PER_BYTE;
    loop {
        let entry = codon_tables.entry(reverse_complement_encoded_sequence, byte_index);
        let first_position = byte_index * POSITIONS_PER_BYTE;
        let last_offset = (scanning_start_position - first_position).min(POSITIONS_PER_BYTE - 1);
        if entry != 0 || first_position <= EDGE_POSITION_LIMIT {
            for offset in (0..=last_offset).rev() {
                let position_index = first_position + offset;
                process_reverse_position(
                    reverse_complement_encoded_sequence,
                    nodes,
                    context,
                    masks,
                    training,
                    position_index,
                    ByteCodonTable::is_stop(entry, offset),
                    ByteCodonTable::start_type(entry, offset),
                );
            }
        }
        if byte_index == 0 {
            break;
        }
        byte_index -= 1;
    }
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn process_reverse_position(
    reverse_complement_encoded_sequence: &[u8],
    nodes: &mut Vec<Node>,
    context: &mut StrandProcessingContext,
    masks: &[Mask],
    training: &Training,
    position_index: usize,
    is_stop_codon: bool,
    start_type: Option<CodonType>,
) {
    {
        let reading_frame_index = position_index % READING_FRAMES;

        if is_stop_codon {
            if context.has_start_codon[reading_frame_index] {
                let node = create_reverse_stop_node(
                    context.last_stop_positions[reading_frame_index],
                    position_index as isize,
                    context.sequence_length,
                    reverse_complement_encoded_sequence,
                    training,
                );
                nodes.push(node);
            }

            context.minimum_distances[reading_frame_index] = MINIMUM_GENE_LENGTH;
            context.last_stop_positions[reading_frame_index] = position_index;
            context.has_start_codon[reading_frame_index] = false;
        } else if context.last_stop_positions[reading_frame_index] < context.sequence_length {
            // A start codon that fails the length or mask checks still falls
            // through to the edge-gene test, as in Prodigal.
            let valid_start = start_type.filter(|_| {
                is_valid_reverse_gene(
                    position_index,
                    context.last_stop_positions[reading_frame_index],
                    context.minimum_distances[reading_frame_index],
                    context.sequence_length,
                    context.circular,
                    masks,
                )
            });
            if let Some(codon_type) = valid_start {
                let node = create_reverse_start_node(
                    position_index,
                    codon_type,
                    context.last_stop_positions[reading_frame_index] as isize,
                    context.sequence_length,
                );

                context.has_start_codon[reading_frame_index] = true;
                nodes.push(node);
            } else if is_reverse_edge_gene(
                position_index,
                context.last_stop_positions[reading_frame_index],
                context.sequence_length,
                context.closed || context.circular,
                masks,
            ) {
                let node = create_reverse_edge_node(
                    position_index,
                    context.last_stop_positions[reading_frame_index] as isize,
                    context.sequence_length,
                );

                context.has_start_codon[reading_frame_index] = true;
                nodes.push(node);
            }
        }
    }
}

/// Handle remaining starts at the end of strand processing
fn handle_remaining_starts(
    encoded_sequence: &[u8],
    nodes: &mut Vec<Node>,
    context: &StrandProcessingContext,
    strand: Strand,
    training: &Training,
) {
    if context.circular {
        return;
    }

    for i in 0..READING_FRAMES {
        if context.has_start_codon[i % READING_FRAMES] {
            let (position_index, stop_value, is_edge) = match strand {
                Strand::Forward => {
                    let is_edge = !is_stop(
                        encoded_sequence,
                        context.last_stop_positions[i % READING_FRAMES],
                        training,
                    );
                    // C code encodes last few forward stops with negative stop_val (i-6)
                    // to allow earliest start connections across sequence edge.
                    let stop_val = (i as isize) - 6;
                    (
                        context.last_stop_positions[i % READING_FRAMES],
                        stop_val,
                        is_edge,
                    )
                }
                Strand::Reverse => {
                    let is_edge = !is_stop(
                        encoded_sequence,
                        context.last_stop_positions[i % READING_FRAMES],
                        training,
                    );
                    let position_index = context.sequence_length
                        - context.last_stop_positions[i % READING_FRAMES]
                        - 1;
                    let stop_val = (context.sequence_length + 5 - i) as isize;
                    (position_index, stop_val, is_edge)
                }
                Strand::Unknown => unreachable!("Unknown strand should not be processed"),
            };

            nodes.push(Node {
                position: NodePosition {
                    index: position_index,
                    codon_type: CodonType::Stop,
                    strand,
                    is_edge,
                    stop_value,
                },
                ..Node::default()
            });
        }
    }
}

/// Create a stop node for forward strand
fn create_stop_node(
    index: usize,
    stop_value: isize,
    strand: Strand,
    encoded_sequence: &[u8],
    training: &Training,
) -> Node {
    let is_edge = !is_stop(encoded_sequence, index, training);

    Node {
        position: NodePosition {
            index,
            codon_type: CodonType::Stop,
            strand,
            is_edge,
            stop_value,
        },
        ..Node::default()
    }
}

/// Create a stop node for reverse strand
fn create_reverse_stop_node(
    index: usize,
    stop_value: isize,
    sequence_length: usize,
    reverse_complement_encoded_sequence: &[u8],
    training: &Training,
) -> Node {
    let is_edge = !is_stop(reverse_complement_encoded_sequence, index, training);

    Node {
        position: NodePosition {
            index: sequence_length - index - 1,
            codon_type: CodonType::Stop,
            strand: Strand::Reverse,
            is_edge,
            stop_value: sequence_length as isize - stop_value - 1,
        },
        ..Node::default()
    }
}

/// Create a start node
fn create_start_node(
    position_index: usize,
    codon_type: CodonType,
    stop_value: isize,
    strand: Strand,
) -> Node {
    Node {
        position: NodePosition {
            index: position_index,
            codon_type,
            strand,
            is_edge: false,
            stop_value,
        },
        ..Node::default()
    }
}

/// Create a reverse start node
fn create_reverse_start_node(
    position: usize,
    codon_type: CodonType,
    stop_value: isize,
    sequence_length: usize,
) -> Node {
    Node {
        position: NodePosition {
            index: sequence_length - position - 1,
            codon_type,
            strand: Strand::Reverse,
            is_edge: false,
            stop_value: sequence_length as isize - stop_value - 1,
        },
        ..Node::default()
    }
}

/// Create an edge node
fn create_edge_node(index: usize, stop_value: isize, strand: Strand) -> Node {
    Node {
        position: NodePosition {
            index,
            codon_type: CodonType::Atg,
            strand,
            is_edge: true,
            stop_value,
        },
        ..Node::default()
    }
}

/// Create a reverse edge node
fn create_reverse_edge_node(position: usize, stop_value: isize, sequence_length: usize) -> Node {
    Node {
        position: NodePosition {
            index: sequence_length - position - 1,
            codon_type: CodonType::Atg,
            strand: Strand::Reverse,
            is_edge: true,
            stop_value: sequence_length as isize - stop_value - 1,
        },
        ..Node::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::encode_sequence;

    fn get_encoded_sequence(input: &[u8]) -> Vec<u8> {
        let sequence_length = input.len();
        let mut seq = vec![0u8; (sequence_length * 2).div_ceil(8)];
        let mut unknown_sequence = vec![0u8; sequence_length.div_ceil(8)];
        let mut masks = Vec::new();
        let _ = encode_sequence(input, &mut seq, &mut unknown_sequence, &mut masks, false).unwrap();
        seq
    }

    fn create_test_encoded_sequence(input: &[u8]) -> EncodedSequence {
        EncodedSequence::without_masking(input)
    }

    fn create_test_training() -> Training {
        Training {
            gc_content: 0.5,
            translation_table: 11,
            uses_shine_dalgarno: true,
            start_type_weights: [1.0, 2.0, 3.0],
            rbs_weights: Box::new([1.0; 28]),
            upstream_composition: Box::new([[0.25; 4]; 32]),
            motif_weights: Box::new([[[1.0; 4096]; 4]; 4]),
            no_motif_weight: 0.5,
            start_weight_factor: 4.35,
            gc_bias_factors: [1.0; 3],
            gene_dicodon_table: Box::new([0.0; 4096]),
            total_dicodons: 0,
        }
    }

    #[test]
    fn test_add_nodes_closed_sequence() {
        let sequence = b"ATGAAATAAGTGAAATAG";
        let encoded_sequence = create_test_encoded_sequence(sequence);
        let mut nodes = Vec::new();
        let training = create_test_training();

        let result = add_nodes(&encoded_sequence, &mut nodes, true, false, &training);

        assert!(result.is_ok());
    }

    #[test]
    fn test_add_nodes_with_masks() {
        let sequence = b"ATGAAATAAGTGAAATAG";
        let encoded_sequence = create_test_encoded_sequence(sequence);
        let mut nodes = Vec::new();
        let training = create_test_training();

        let result = add_nodes(&encoded_sequence, &mut nodes, false, false, &training);

        assert!(result.is_ok());
    }

    #[test]
    fn test_initialize_strand_arrays() {
        let sequence_length = 100;
        let sequence_length_mod = 1;
        let closed = false;
        let circular = false;

        let context =
            StrandProcessingContext::new(sequence_length, sequence_length_mod, closed, circular);

        // Check that arrays were properly initialized
        assert_eq!(context.has_start_codon, [false; READING_FRAMES]);
        for &dist in &context.minimum_distances {
            assert_eq!(dist, MINIMUM_EDGE_GENE_LENGTH);
        }
        // Let's just verify they are reasonable values
        for &val in &context.last_stop_positions {
            assert!(val <= sequence_length + READING_FRAMES);
        }
    }

    #[test]
    fn test_initialize_strand_arrays_closed() {
        let sequence_length = 10;
        let sequence_length_mod = 0;
        let closed = true;
        let circular = false;

        let context =
            StrandProcessingContext::new(sequence_length, sequence_length_mod, closed, circular);

        // Check that arrays were properly initialized
        assert_eq!(context.has_start_codon, [false; READING_FRAMES]);
        for &dist in &context.minimum_distances {
            assert_eq!(dist, MINIMUM_EDGE_GENE_LENGTH);
        }
    }

    #[test]
    fn test_create_start_node() {
        let node = create_start_node(100, CodonType::Atg, 200, Strand::Forward);

        assert_eq!(node.position.index, 100);
        assert_eq!(node.position.codon_type, CodonType::Atg);
        assert_eq!(node.position.strand, Strand::Forward);
        assert_eq!(node.position.stop_value, 200);
        assert!(!node.position.is_edge);
    }

    #[test]
    fn test_create_reverse_start_node() {
        let sequence_length = 1000;
        let node = create_reverse_start_node(100, CodonType::Gtg, 200, sequence_length);

        assert_eq!(node.position.index, sequence_length - 100 - 1);
        assert_eq!(node.position.codon_type, CodonType::Gtg);
        assert_eq!(node.position.strand, Strand::Reverse);
        assert_eq!(
            node.position.stop_value,
            (sequence_length - 200 - 1) as isize
        );
        assert!(!node.position.is_edge);
    }

    #[test]
    fn test_create_edge_node() {
        let node = create_edge_node(50, 150, Strand::Forward);

        assert_eq!(node.position.index, 50);
        assert_eq!(node.position.codon_type, CodonType::Atg);
        assert_eq!(node.position.strand, Strand::Forward);
        assert_eq!(node.position.stop_value, 150);
        assert!(node.position.is_edge);
    }

    #[test]
    fn test_create_reverse_edge_node() {
        let sequence_length = 1000;
        let node = create_reverse_edge_node(50, 150, sequence_length);

        assert_eq!(node.position.index, sequence_length - 50 - 1);
        assert_eq!(node.position.codon_type, CodonType::Atg);
        assert_eq!(node.position.strand, Strand::Reverse);
        assert_eq!(
            node.position.stop_value,
            (sequence_length - 150 - 1) as isize
        );
        assert!(node.position.is_edge);
    }

    #[test]
    fn test_create_stop_node() {
        let sequence = b"TAAGGG";
        let encoded_seq = get_encoded_sequence(sequence);
        let training = create_test_training();

        let node = create_stop_node(0, 3, Strand::Forward, &encoded_seq, &training);

        assert_eq!(node.position.index, 0);
        assert_eq!(node.position.codon_type, CodonType::Stop);
        assert_eq!(node.position.strand, Strand::Forward);
        assert_eq!(node.position.stop_value, 3);
    }

    #[test]
    fn test_create_reverse_stop_node() {
        let sequence = b"TAAGGG";
        let encoded_seq = get_encoded_sequence(sequence);
        let training = create_test_training();
        let sequence_length = sequence.len();

        let node = create_reverse_stop_node(0, 3, sequence_length, &encoded_seq, &training);

        assert_eq!(node.position.index, sequence_length - 1);
        assert_eq!(node.position.codon_type, CodonType::Stop);
        assert_eq!(node.position.strand, Strand::Reverse);
        assert_eq!(node.position.stop_value, (sequence_length - 3 - 1) as isize);
    }

    #[test]
    fn test_handle_remaining_starts_forward() {
        let sequence = b"ATGAAA";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut nodes = Vec::new();
        let training = create_test_training();

        // Create a context with some start codons detected
        let mut context = StrandProcessingContext::new(sequence.len(), 0, false, false);
        context.last_stop_positions = [0, READING_FRAMES, 6];
        context.has_start_codon = [true, false, true];

        handle_remaining_starts(
            &encoded_seq,
            &mut nodes,
            &context,
            Strand::Forward,
            &training,
        );

        assert_eq!(nodes.len(), 2);
        assert_eq!(nodes[0].position.strand, Strand::Forward);
        assert_eq!(nodes[1].position.strand, Strand::Forward);
    }

    #[test]
    fn test_handle_remaining_starts_reverse() {
        let sequence = b"ATGAAA";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut nodes = Vec::new();
        let training = create_test_training();

        // Create a context with some start codons detected
        let mut context = StrandProcessingContext::new(sequence.len(), 0, false, false);
        context.last_stop_positions = [0, READING_FRAMES, 6];
        context.has_start_codon = [false, true, false];

        handle_remaining_starts(
            &encoded_seq,
            &mut nodes,
            &context,
            Strand::Reverse,
            &training,
        );

        assert_eq!(nodes.len(), 1);
        assert_eq!(nodes[0].position.strand, Strand::Reverse);
    }

    #[test]
    fn test_handle_remaining_starts_no_starts() {
        let sequence = b"ATGAAA";
        let encoded_seq = get_encoded_sequence(sequence);
        let mut nodes = Vec::new();
        let training = create_test_training();

        // Create a context with no start codons detected
        let mut context = StrandProcessingContext::new(sequence.len(), 0, false, false);
        context.last_stop_positions = [0, READING_FRAMES, 6];
        context.has_start_codon = [false, false, false];

        handle_remaining_starts(
            &encoded_seq,
            &mut nodes,
            &context,
            Strand::Forward,
            &training,
        );

        assert!(nodes.is_empty());
    }

    #[test]
    fn test_add_nodes_short_sequence() {
        let sequence = b"ATG";
        let encoded_sequence = create_test_encoded_sequence(sequence);
        let mut nodes = Vec::new();
        let training = create_test_training();

        let result = add_nodes(&encoded_sequence, &mut nodes, false, false, &training);

        assert!(result.is_ok());
    }

    #[test]
    fn test_add_nodes_circular_detects_wrapped_forward_start() {
        // Build a sequence with a wrapped forward ORF that clears minimum gene length.
        let mut sequence = vec![b'C'; 150];
        sequence[20..23].copy_from_slice(b"TAA");
        sequence[80..83].copy_from_slice(b"ATG");

        let encoded_sequence = create_test_encoded_sequence(&sequence);
        let training = create_test_training();

        let mut linear_nodes = Vec::new();
        let _ = add_nodes(&encoded_sequence, &mut linear_nodes, true, false, &training).unwrap();

        let mut circular_nodes = Vec::new();
        let _ = add_nodes(
            &encoded_sequence,
            &mut circular_nodes,
            false,
            true,
            &training,
        )
        .unwrap();

        let has_wrapped_forward_start_linear = linear_nodes.iter().any(|n| {
            n.position.strand == Strand::Forward
                && n.position.codon_type != CodonType::Stop
                && n.position.index == 80
                && n.position.stop_value == 20
        });
        let has_wrapped_forward_start_circular = circular_nodes.iter().any(|n| {
            n.position.strand == Strand::Forward
                && n.position.codon_type != CodonType::Stop
                && n.position.index == 80
                && n.position.stop_value == 20
        });

        assert!(!has_wrapped_forward_start_linear);
        assert!(has_wrapped_forward_start_circular);
    }
}
