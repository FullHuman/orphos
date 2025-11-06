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
}

impl StrandProcessingContext {
    /// Create a new processing context with initialized arrays
    fn new(sequence_length: usize, sequence_length_mod: usize, closed: bool) -> Self {
        let mut context = Self {
            last_stop_positions: ReadingFrameArray::default(),
            has_start_codon: ReadingFrameArray::default(),
            minimum_distances: ReadingFrameArray::default(),
            sequence_length,
            closed,
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
}

/// Add nodes for start and stop codons in both directions
///
/// This function identifies potential start and stop codons in the sequence
/// and creates corresponding nodes for gene prediction analysis.
///
/// # Arguments
/// * `encoded_sequence` - The complete encoded sequence data
/// * `nodes` - Vector to store the created nodes
/// * `closed` - Whether to treat sequence as circular/closed
/// * `training` - Training data for scoring parameters
///
/// # Returns
/// Number of nodes created, or error if processing fails
pub fn add_nodes(
    encoded_sequence: &EncodedSequence,
    nodes: &mut Vec<Node>,
    closed: bool,
    training: &Training,
) -> Result<usize, OrphosError> {
    let sequence_length = encoded_sequence.sequence_length;
    // Clear the nodes vector
    nodes.clear();
    nodes.reserve(STT_NOD);

    let sequence_length_mod = sequence_length % READING_FRAMES;

    let mut context = StrandProcessingContext::new(sequence_length, sequence_length_mod, closed);
    process_forward_strand(
        &encoded_sequence.forward_sequence,
        nodes,
        &mut context,
        &encoded_sequence.masks,
        training,
    );
    handle_remaining_starts(
        &encoded_sequence.forward_sequence,
        nodes,
        &context,
        Strand::Forward,
        training,
    );

    let mut context = StrandProcessingContext::new(sequence_length, sequence_length_mod, closed);
    process_reverse_strand(
        &encoded_sequence.reverse_complement_sequence,
        nodes,
        &mut context,
        &encoded_sequence.masks,
        training,
    );
    handle_remaining_starts(
        &encoded_sequence.reverse_complement_sequence,
        nodes,
        &context,
        Strand::Reverse,
        training,
    );

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
) {
    let scanning_start_position = context.sequence_length.saturating_sub(CODON_LENGTH);

    for position_index in (0..=scanning_start_position).rev() {
        let reading_frame_index = position_index % READING_FRAMES;

        if is_stop(encoded_sequence, position_index, training) {
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
            continue;
        }

        if context.last_stop_positions[reading_frame_index] >= context.sequence_length {
            continue;
        }

        // Check for start codons with unified logic
        if let Some(codon_type) = check_start_codon(encoded_sequence, position_index, training) {
            if is_valid_gene(
                position_index,
                context.last_stop_positions[reading_frame_index],
                context.minimum_distances[reading_frame_index],
                masks,
            ) {
                let node = create_start_node(
                    position_index,
                    codon_type,
                    context.last_stop_positions[reading_frame_index] as isize,
                    Strand::Forward,
                );
                context.has_start_codon[reading_frame_index] = true;
                nodes.push(node);
            }
        } else if is_edge_gene(
            position_index,
            context.last_stop_positions[reading_frame_index],
            context.closed,
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

/// Process reverse strand to find genes
fn process_reverse_strand(
    reverse_complement_encoded_sequence: &[u8],
    nodes: &mut Vec<Node>,
    context: &mut StrandProcessingContext,
    masks: &[Mask],
    training: &Training,
) {
    let scanning_start_position = context.sequence_length.saturating_sub(CODON_LENGTH);

    for position_index in (0..=scanning_start_position).rev() {
        let reading_frame_index = position_index % READING_FRAMES;

        if is_stop(
            reverse_complement_encoded_sequence,
            position_index,
            training,
        ) {
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
            continue;
        }

        if context.last_stop_positions[reading_frame_index] >= context.sequence_length {
            continue;
        }

        // Check for start codons on reverse strand
        if let Some(codon_type) = check_start_codon(
            reverse_complement_encoded_sequence,
            position_index,
            training,
        ) {
            if is_valid_reverse_gene(
                position_index,
                context.last_stop_positions[reading_frame_index],
                context.minimum_distances[reading_frame_index],
                context.sequence_length,
                masks,
            ) {
                let node = create_reverse_start_node(
                    position_index,
                    codon_type,
                    context.last_stop_positions[reading_frame_index] as isize,
                    context.sequence_length,
                );

                context.has_start_codon[reading_frame_index] = true;
                nodes.push(node);
            }
        } else if is_reverse_edge_gene(
            position_index,
            context.last_stop_positions[reading_frame_index],
            context.sequence_length,
            context.closed,
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

/// Handle remaining starts at the end of strand processing
fn handle_remaining_starts(
    encoded_sequence: &[u8],
    nodes: &mut Vec<Node>,
    context: &StrandProcessingContext,
    strand: Strand,
    training: &Training,
) {
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

        let result = add_nodes(&encoded_sequence, &mut nodes, true, &training);

        assert!(result.is_ok());
    }

    #[test]
    fn test_add_nodes_with_masks() {
        let sequence = b"ATGAAATAAGTGAAATAG";
        let encoded_sequence = create_test_encoded_sequence(sequence);
        let mut nodes = Vec::new();
        let training = create_test_training();

        let result = add_nodes(&encoded_sequence, &mut nodes, false, &training);

        assert!(result.is_ok());
    }

    #[test]
    fn test_initialize_strand_arrays() {
        let sequence_length = 100;
        let sequence_length_mod = 1;
        let closed = false;

        let context = StrandProcessingContext::new(sequence_length, sequence_length_mod, closed);

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

        let context = StrandProcessingContext::new(sequence_length, sequence_length_mod, closed);

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
        let mut context = StrandProcessingContext::new(sequence.len(), 0, false);
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
        let mut context = StrandProcessingContext::new(sequence.len(), 0, false);
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
        let mut context = StrandProcessingContext::new(sequence.len(), 0, false);
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

        let result = add_nodes(&encoded_sequence, &mut nodes, false, &training);

        assert!(result.is_ok());
    }
}
