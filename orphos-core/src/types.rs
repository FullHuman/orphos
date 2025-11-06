use std::fmt;

use bio::bio_types::strand::Strand;
use thiserror::Error;

/// Weights for start codon recognition (ATG, GTG, TTG).
///
/// Array of three weights indicating the relative preference for each
/// start codon type based on training data.
pub type StartWeights = [f64; 3];

/// Weights for ribosome binding site scoring (28 possible motifs).
///
/// Scores for different RBS patterns discovered during training.
pub type RbsWeights = [f64; 28];

/// Upstream sequence composition frequencies \[position\]\[nucleotide\].
///
/// Nucleotide frequencies at positions upstream of start codons,
/// used for start site scoring.
pub type UpstreamComposition = [[f64; 4]; 32];

/// Motif scoring weights \[motif_type\]\[len\]\[pattern\].
///
/// Weights for non-SD start recognition motifs.
pub type MotifWeights = [[[f64; 4096]; 4]; 4];

/// GC content bias scores for each reading frame.
///
/// Adjustment factors for each frame based on GC content analysis.
pub type GcFrameScores = [f64; 3];

/// Pointers to best scoring nodes in each frame.
///
/// Used in dynamic programming to track optimal paths.
pub type StarPointers = [Option<usize>; 3];

/// Dicodon frequency table for coding potential scoring.
///
/// 4096 entries (64x64) for all possible dicodon pairs, used to
/// calculate coding vs non-coding likelihood.
pub type DicodonTable = [f64; 4096];

/// GC bias adjustment factors for each reading frame.
///
/// Multipliers to correct for GC-biased reading frames.
pub type GcBias = [f64; 3];

/// Types of start codons recognized by Orphos.
///
/// Prokaryotic genes primarily use three start codons: ATG (most common),
/// GTG, and TTG. Stop codons are also tracked for ORF detection.
///
/// # Examples
///
/// ```rust
/// use orphos_core::types::CodonType;
///
/// let atg = CodonType::Atg;
/// assert_eq!(atg.to_index(), 0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CodonType {
    /// ATG start codon
    Atg,
    /// GTG start codon
    Gtg,
    /// TTG start codon
    Ttg,
    /// Stop codon (TAA, TAG, TGA)
    Stop,
}

impl CodonType {
    /// Convert codon type to array index for scoring
    #[must_use]
    pub const fn to_index(self) -> usize {
        match self {
            Self::Atg => 0,
            Self::Gtg => 1,
            Self::Ttg => 2,
            Self::Stop => 3,
        }
    }
}

/// Gene start types including edge cases
#[derive(Debug, Clone, Default)]
pub enum StartType {
    /// ATG start codon
    Atg,
    /// GTG start codon
    Gtg,
    /// TTG start codon
    Ttg,
    /// Gene extends to sequence edge
    Edge,
    /// Unknown or unrecognized start
    #[default]
    Unknown,
}

impl fmt::Display for StartType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Atg => write!(f, "ATG"),
            Self::Gtg => write!(f, "GTG"),
            Self::Ttg => write!(f, "TTG"),
            Self::Edge => write!(f, "Edge"),
            Self::Unknown => write!(f, "Unknown"),
        }
    }
}

impl From<CodonType> for StartType {
    fn from(codon_type: CodonType) -> Self {
        match codon_type {
            CodonType::Atg => Self::Atg,
            CodonType::Gtg => Self::Gtg,
            CodonType::Ttg => Self::Ttg,
            CodonType::Stop => Self::Unknown,
        }
    }
}

/// Ribosome binding site motif information
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Motif {
    /// Index of the motif in the motif table
    pub index: usize,
    /// Length of the motif sequence
    pub length: usize,
    /// Index indicating spacer position
    pub space_index: usize,
    /// Distance between motif and start codon
    pub spacer: usize,
    /// Calculated motif score
    pub score: f64,
}

/// Positional information for a gene prediction node
#[derive(Debug, Clone)]
pub struct NodePosition {
    /// Position index in the sequence
    pub index: usize,
    /// Strand orientation (forward or reverse)
    pub strand: Strand,
    /// Type of codon at this position
    pub codon_type: CodonType,
    /// Position of the corresponding stop codon
    pub stop_value: isize,
    /// Whether this gene extends to the sequence edge
    pub is_edge: bool,
}

impl Default for NodePosition {
    fn default() -> Self {
        Self {
            index: 0,
            strand: Strand::Forward,
            codon_type: CodonType::Atg,
            stop_value: 0,
            is_edge: false,
        }
    }
}

/// Gene scoring information for dynamic programming
#[derive(Debug, Clone, Default)]
pub struct NodeScores {
    /// GC content of the predicted gene
    pub gc_content: f64,
    /// Coding potential score based on dicodon frequencies
    pub coding_score: f64,
    /// Start codon recognition score
    pub start_score: f64,
    /// Ribosome binding site score
    pub ribosome_binding_score: f64,
    /// Start codon type preference score (ATG vs GTG vs TTG)
    pub type_score: f64,
    /// Upstream sequence composition score
    pub upstream_score: f64,
    /// Total combined score for gene calling
    pub total_score: f64,
    /// GC content bias scores for each reading frame
    pub gc_frame_scores: GcFrameScores,
}

/// Dynamic programming state information for gene finding
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct NodeState {
    /// Pointers to optimal start nodes in each reading frame
    pub start_pointers: StarPointers,
    /// Traceback pointer for dynamic programming path reconstruction
    pub traceback: Option<usize>,
    /// Forward trace pointer for gene path tracking
    pub trace_forward: Option<usize>,
    /// Marker for overlapping gene resolution
    pub overlap_marker: Option<usize>,
    /// Whether this node was eliminated during gene selection
    pub is_eliminated: bool,
    /// Reading frame with highest GC bias score
    pub gc_bias_frame: usize,
}

/// Ribosome binding site and motif information for a node
#[derive(Debug, Clone, Default)]
pub struct NodeMotifInfo {
    /// Indices of ribosome binding site motifs [exact, mismatch]
    pub ribosome_binding_sites: [usize; 2],
    /// Best scoring upstream motif for start recognition
    pub best_motif: Motif,
}

/// A gene prediction node containing position, scoring, and state information
#[derive(Debug, Clone, Default)]
pub struct Node {
    /// Physical position and codon information
    pub position: NodePosition,
    /// Calculated scores for gene quality assessment
    pub scores: NodeScores,
    /// Dynamic programming state for optimal path finding
    pub state: NodeState,
    /// Ribosome binding site and motif information
    pub motif_info: NodeMotifInfo,
}

pub trait NodeSliceExtension {
    /// Iterate through the traceback path starting from a given index
    fn iter_traceback(&self, start_index: usize) -> impl Iterator<Item = &Node>;

    /// Iterate through the forward trace path starting from a given index
    fn iter_trace_forward(&self, start_index: usize) -> impl Iterator<Item = &Node>;

    /// Iterate through the traceback path starting from a given index, returning indices
    fn iter_traceback_indices(&self, start_index: usize) -> impl Iterator<Item = usize>;

    /// Iterate through the forward trace path starting from a given index, returning indices
    fn iter_trace_forward_indices(&self, start_index: usize) -> impl Iterator<Item = usize>;
}

impl NodeSliceExtension for [Node] {
    fn iter_traceback(&self, start_index: usize) -> impl Iterator<Item = &Node> {
        std::iter::successors(self.get(start_index), |previous_node| {
            previous_node
                .state
                .traceback
                .and_then(|path_index| self.get(path_index))
        })
    }

    fn iter_trace_forward(&self, start_index: usize) -> impl Iterator<Item = &Node> {
        std::iter::successors(self.get(start_index), |previous_node| {
            previous_node
                .state
                .trace_forward
                .and_then(|path_index| self.get(path_index))
        })
    }

    fn iter_traceback_indices(&self, start_index: usize) -> impl Iterator<Item = usize> {
        std::iter::successors(Some(start_index), |&previous_index| {
            self.get(previous_index)
                .and_then(|node| node.state.traceback)
        })
    }

    fn iter_trace_forward_indices(&self, start_index: usize) -> impl Iterator<Item = usize> {
        std::iter::successors(Some(start_index), |&previous_index| {
            self.get(previous_index)
                .and_then(|node| node.state.trace_forward)
        })
    }
}

/// Training data structure containing all learned parameters
#[derive(Debug)]
pub struct Training {
    /// GC content of the sequence
    pub gc_content: f64,
    /// Translation table (11 for bacteria, 4 for archaea, etc.)
    pub translation_table: i32,
    /// Whether organism uses Shine-Dalgarno motifs
    pub uses_shine_dalgarno: bool,
    /// Start codon type weights [ATG, GTG, TTG]
    pub start_type_weights: StartWeights,
    /// RBS/Shine-Dalgarno motif weights (28 different patterns)
    pub rbs_weights: Box<RbsWeights>,
    /// Upstream composition weights \[position\]\[base\]
    pub upstream_composition: Box<UpstreamComposition>,
    /// Non-SD motif weights \[length-3\]\[spacer\]\[motif_index\]
    pub motif_weights: Box<MotifWeights>,
    /// Weight for no motif found
    pub no_motif_weight: f64,
    /// Start weight factor
    pub start_weight_factor: f64,
    /// GC bias factors
    pub gc_bias_factors: GcBias,
    /// Dicodon gene scoring table
    pub gene_dicodon_table: Box<DicodonTable>,
    /// Total dicodons counted across genes (for parity/debugging)
    pub total_dicodons: u32,
}

impl Default for Training {
    fn default() -> Self {
        Self {
            gc_content: 0.5,
            translation_table: 11,
            uses_shine_dalgarno: true,
            start_type_weights: [0.0; 3],
            rbs_weights: Box::new([0.0; 28]),
            upstream_composition: Box::new([[0.0; 4]; 32]),
            motif_weights: {
                let vec = vec![[[0.0; 4096]; 4]; 4];
                vec.into_boxed_slice().try_into().unwrap()
            },
            no_motif_weight: 0.0,
            start_weight_factor: 4.35,
            gc_bias_factors: [1.0; 3],
            gene_dicodon_table: Box::new([0.0; 4096]),
            total_dicodons: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GeneCoordinates {
    /// Start position of the gene (1-based)
    pub begin: usize,
    /// End position of the gene (1-based, inclusive)
    pub end: usize,
    /// Strand orientation (forward/reverse)
    pub strand: Strand,
    /// Sequence index of the start codon
    pub start_index: usize,
    /// Sequence index of the stop codon
    pub stop_index: usize,
}

impl Default for GeneCoordinates {
    fn default() -> Self {
        Self {
            begin: 0,
            end: 0,
            strand: Strand::Unknown,
            start_index: 0,
            stop_index: 0,
        }
    }
}

/// Gene annotation metadata including identification and analysis details
#[derive(Debug, Clone, Default)]
pub struct GeneAnnotation {
    /// Unique identifier for the gene
    pub identifier: String,
    /// Whether the gene is truncated at the 5' end (partial left)
    pub is_partial_left: bool,
    /// Whether the gene is truncated at the 3' end (partial right)
    pub is_partial_right: bool,
    /// Type of start codon used
    pub start_type: StartType,
    /// Ribosome binding site motif sequence, if found
    pub ribosome_binding_motif: Option<String>,
    /// Spacer sequence between RBS and start codon, if applicable
    pub ribosome_binding_spacer: Option<String>,
    /// GC content of the predicted gene
    pub gc_content: f64,
}
impl GeneAnnotation {
    /// Create a new gene annotation with basic information
    #[must_use]
    pub const fn new(
        id: String,
        partial_left: bool,
        partial_right: bool,
        start_type: StartType,
        gc_content: f64,
    ) -> Self {
        Self {
            identifier: id,
            is_partial_left: partial_left,
            is_partial_right: partial_right,
            start_type,
            ribosome_binding_motif: None,
            ribosome_binding_spacer: None,
            gc_content,
        }
    }

    /// Add ribosome binding site information to the annotation
    #[must_use]
    pub fn with_rbs(mut self, motif: String, spacer: String) -> Self {
        self.ribosome_binding_motif = Some(motif);
        self.ribosome_binding_spacer = Some(spacer);
        self
    }

    /// Get RBS motif for display, returning "None" if not present
    fn rbs_motif_display(&self) -> &str {
        self.ribosome_binding_motif.as_deref().unwrap_or("None")
    }

    /// Get RBS spacer for display, returning "None" if not present
    fn rbs_spacer_display(&self) -> &str {
        self.ribosome_binding_spacer.as_deref().unwrap_or("None")
    }
}
impl fmt::Display for GeneAnnotation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "ID={};partial={}{};start_type={};rbs_motif={};rbs_spacer={};gc_cont={:.3}",
            self.identifier,
            if self.is_partial_left { "1" } else { "0" },
            if self.is_partial_right { "1" } else { "0" },
            self.start_type,
            self.rbs_motif_display(),
            self.rbs_spacer_display(),
            self.gc_content
        )
    }
}

/// Comprehensive scoring information for a predicted gene
#[derive(Debug, Clone, Default)]
pub struct GeneScore {
    /// Overall confidence score for the gene prediction
    pub confidence: f64,
    /// Total combined score from all components
    pub total_score: f64,
    /// Coding potential score based on sequence composition
    pub coding_score: f64,
    /// Start codon recognition score
    pub start_score: f64,
    /// Ribosome binding site score
    pub ribosome_binding_score: f64,
    /// Upstream sequence composition score
    pub upstream_score: f64,
    /// Start codon type preference score (ATG vs GTG vs TTG)
    pub type_score: f64,
}

impl fmt::Display for GeneScore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "conf={:.2};score={:.2};cscore={:.2};sscore={:.2};rscore={:.2};uscore={:.2};tscore={:.2}",
            self.confidence,
            self.total_score,
            self.coding_score,
            self.start_score,
            self.ribosome_binding_score,
            self.upstream_score,
            self.type_score
        )
    }
}

/// Complete gene prediction result with coordinates, scores, and annotations
#[derive(Debug, Clone, Default)]
pub struct Gene {
    /// Physical coordinates and strand information
    pub coordinates: GeneCoordinates,
    /// Detailed scoring information for quality assessment
    pub score: GeneScore,
    /// Metadata and annotation details
    pub annotation: GeneAnnotation,
}

/// Mask for regions to exclude from gene prediction (e.g., runs of N characters)
#[derive(Debug, Clone)]
pub struct Mask {
    /// Start position of the masked region
    pub begin: usize,
    /// End position of the masked region
    pub end: usize,
}

/// Error types that can occur during gene prediction analysis
#[derive(Error, Debug)]
pub enum OrphosError {
    /// Invalid input sequence format or content
    #[error("Invalid sequence: {0}")]
    InvalidSequence(String),
    /// Problem with training file format or data
    #[error("Invalid training file: {0}")]
    InvalidTrainingFile(String),
    /// File I/O operation failed
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    /// Error parsing input data
    #[error("Parse error: {0}")]
    ParseError(String),
    /// Invalid strand specification
    #[error("Invalid strand")]
    InvalidStrand,
    /// Sequence length is invalid for analysis
    #[error("Invalid sequence length")]
    InvalidSequenceLength,
    /// Node position is out of sequence bounds
    #[error("Invalid node position: {0} (sequence length: {1})")]
    InvalidNodePosition(usize, usize),
    /// Motif position is invalid
    #[error("Invalid motif position")]
    InvalidMotifPosition,
    /// Invalid spacer distance
    #[error("Invalid spacer")]
    InvalidSpacer,
    /// Motif index is out of range
    #[error("Invalid motif index")]
    InvalidMotifIndex,
}
