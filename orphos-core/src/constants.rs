// =============================================================================
// =============================================================================

// =============================================================================
// =============================================================================

/// Version string for Orphos
pub const VERSION: &str = "1.0.0";

// =============================================================================
// =============================================================================

/// Maximum allowed sequence length in base pairs
pub const MAX_SEQUENCE_LENGTH: usize = 32_000_000;

/// Maximum line length for input parsing
pub const MAX_LINE_LENGTH: usize = 10_000;

/// Size of sequence masks for filtering regions
pub const MASK_SIZE: usize = 50;

/// Minimum sequence length required for gene prediction
pub const MIN_SEQUENCE_LENGTH: usize = 96;

/// Number of reading frames to analyze (forward and reverse)
pub const READING_FRAMES: usize = 3;

/// Length of a codon in base pairs
pub const CODON_LENGTH: usize = 3;

/// Sliding window size for sequence analysis
pub const SLIDING_WINDOW_SIZE: usize = 120;

/// SIMD processing chunk size for sequence analysis
pub const CHUNK_SIZE: usize = 32;

// =============================================================================
// =============================================================================

/// Minimum gene length in base pairs for normal genes
pub const MINIMUM_GENE_LENGTH: usize = 90;

/// Minimum gene length in base pairs for genes extending to sequence edges
pub const MINIMUM_EDGE_GENE_LENGTH: usize = 60;

/// Maximum allowed overlap between genes on the same strand
pub const MAXIMUM_SAME_OVERLAP: usize = 60;

/// Maximum allowed overlap between genes on opposite strands
pub const MAXIMUM_OPPOSITE_OVERLAP: i32 = 200;

/// Maximum overlap for optimization algorithms
pub const MAXIMUM_OVERLAP: usize = 60;

/// Offset for stop codon processing
pub const STOP_CODON_OFFSET: usize = 2;

/// Penalty factor for gene overlaps
pub const OVERLAP_PENALTY_FACTOR: f64 = 0.15;

// =============================================================================
// =============================================================================

/// Maximum number of nodes (start/stop positions) to consider
pub const STT_NOD: usize = 100_000;

/// Maximum number of genes that can be predicted in a single sequence
pub const MAXIMUM_GENES: usize = 50000;

/// Alternative maximum gene count constant (used in different contexts)
pub const MAX_GENES: usize = 30000;

/// Maximum distance between nodes for scoring
pub const MAX_NODE_DIST: usize = 500;

// =============================================================================
// =============================================================================

/// Size of dicodon sequences for gene scoring
pub const DICODON_SIZE: usize = 6;

/// Total number of possible dicodons (4^6)
pub const NUM_DICODONS: usize = 4096;

/// Minimum motif length for pattern matching
pub const MIN_MOTIF_LENGTH: usize = 3;

/// Maximum motif length for pattern matching
pub const MAX_MOTIF_LENGTH: usize = 6;

/// Minimum distance from gene start for motif searching
pub const MIN_DISTANCE_FROM_START: usize = 4;

/// Maximum distance for ribosome binding site detection
pub const MAX_RIBOSOME_DISTANCE: usize = 15;

/// Upstream distance for RBS motif searching
pub const RBS_UPSTREAM_DISTANCE: usize = 20;

/// Downstream distance for RBS motif searching
pub const RBS_DOWNSTREAM_DISTANCE: usize = 6;

/// Distance threshold for operon structure consideration
pub const OPERON_DISTANCE: f64 = 60.0;

// =============================================================================
// =============================================================================

/// Minimum log likelihood value for scoring
pub const MIN_LOG_LIKELIHOOD: f64 = -5.0;

/// Maximum log likelihood value for scoring
pub const MAX_LOG_LIKELIHOOD: f64 = 5.0;

/// Minimum motif score threshold
pub const MIN_MOTIF_SCORE: f64 = -4.0;

/// Initial maximum score for optimization
pub const INITIAL_MAX_SCORE: f64 = -100.0;

/// Motif threshold offset for scoring adjustments
pub const MOTIF_THRESHOLD_OFFSET: f64 = 0.69;

/// Minimum cumulative score for sequence analysis
pub const MIN_CUMULATIVE_SCORE: f64 = 6.0;

/// Score bonus applied to genes extending to sequence edges
pub const EDGE_BONUS: f64 = 0.74;

/// Upstream score penalty for edge genes
pub const EDGE_UPSTREAM: f64 = -1.00;

/// Score penalty applied in metagenomic mode
pub const METAGENOMIC_PENALTY: f64 = 7.5;

/// Coefficient divisor for metagenomic penalty calculation
pub const METAGENOMIC_PENALTY_DIVISOR: f64 = 2700.0;

/// Sentinel value for coding scores
pub const CODING_SCORE_SENTINEL: f64 = -10000.0;

// =============================================================================
// =============================================================================

/// Search window size for node analysis
pub const NODE_SEARCH_WINDOW: usize = 500;

/// Alternative search window for optimization
pub const SEARCH_WINDOW: usize = 100;

/// Threshold for short gene classification
pub const SHORT_GENE_THRESHOLD: usize = 250;

/// Length threshold for metagenomic analysis
pub const METAGENOMIC_LENGTH_THRESHOLD: usize = 3000;

/// Minimum gene length for metagenomic analysis
pub const MIN_META_GENE_LENGTH: usize = 120;

/// Length factor for metagenomic minimum calculation
pub const METAGENOMIC_MIN_LENGTH_FACTOR: usize = 1500;

/// Range for upstream motif scanning
pub const UPSTREAM_SCAN_RANGE: usize = 45;

/// Starting position to skip in upstream scanning
pub const UPSTREAM_SKIP_START: usize = 2;

/// Ending position to skip in upstream scanning
pub const UPSTREAM_SKIP_END: usize = 15;

/// Threshold for coding score evaluation
pub const CODING_SCORE_THRESHOLD: f64 = 5.0;

/// Length factor threshold for gene scoring
pub const LENGTH_FACTOR_THRESHOLD: f64 = 3.0;

/// Multiplier for length factor calculations
pub const LENGTH_FACTOR_MULTIPLIER: f64 = 0.5;

/// Threshold for edge position detection
pub const EDGE_POSITION_THRESHOLD: usize = 2;

/// Offset for edge position calculations
pub const EDGE_POSITION_OFFSET: usize = 3;

/// Weight factor for upstream composition scoring
pub const UPSTREAM_COMPOSITION_WEIGHT: f64 = 0.4;

/// Threshold for no motif detection
pub const NO_MOTIF_THRESHOLD: f64 = -0.5;

/// Penalty for negative scores
pub const NEGATIVE_SCORE_PENALTY: f64 = 0.5;

/// Minimum gene size in codons
pub const MIN_GENE_SIZE_CODONS: i32 = 80;

/// Maximum gene size in codons
pub const MAX_GENE_SIZE_CODONS: i32 = 1000;

/// Scaling factor for gene size calculations
pub const GENE_SIZE_SCALING_FACTOR: f64 = 920.0;

// =============================================================================
// =============================================================================

/// Maximum training iterations for Shine-Dalgarno analysis
pub const MAX_TRAINING_ITERATIONS_SD: usize = 10;

/// Maximum training iterations for non-Shine-Dalgarno analysis
pub const MAX_TRAINING_ITERATIONS_NONSD: usize = 20;

/// Initial score threshold for training
pub const INITIAL_SCORE_THRESHOLD: f64 = 35.0;

/// Threshold divisor for training iterations
pub const THRESHOLD_DIVISOR: f64 = 2.0;

/// Gene ratio threshold for training validation
pub const GENE_RATIO_THRESHOLD: f64 = 2000.0;

/// Coverage threshold for upstream motif analysis
pub const UPSTREAM_MOTIF_COVERAGE_THRESHOLD: f64 = 0.2;

/// Minimum GC content for analysis
pub const MIN_GC_CONTENT: f64 = 0.1;

/// Maximum GC content for analysis
pub const MAX_GC_CONTENT: f64 = 0.9;

/// Low GC frequency threshold
pub const LOW_GC_FREQ: f64 = 0.45;

/// High GC frequency threshold
pub const HIGH_GC_FREQ: f64 = 0.05;

/// Minimum weight clamping value
pub const WEIGHT_CLAMP_MIN: f64 = -4.0;

/// Maximum weight clamping value
pub const WEIGHT_CLAMP_MAX: f64 = 4.0;

/// High threshold for RBS weight evaluation
pub const RBS_WEIGHT_THRESHOLD_HIGH: f64 = 1.0;

/// Low threshold for RBS weight evaluation
pub const RBS_WEIGHT_THRESHOLD_LOW: f64 = -0.5;

/// Strong threshold for RBS weight evaluation
pub const RBS_WEIGHT_STRONG_THRESHOLD: f64 = 2.0;

/// Number of codon types (ATG, GTG, TTG)
pub const NUM_CODON_TYPES: usize = 3;

/// Number of RBS weight patterns
pub const NUM_RBS_WEIGHTS: usize = 28;

/// Number of nucleotide bases (A, T, G, C)
pub const NUM_BASES: usize = 4;

/// Number of upstream positions for analysis
pub const UPSTREAM_POSITIONS: usize = 32;

/// Number of motif size categories
pub const NUM_MOTIF_SIZES: usize = 4;

/// Maximum motif index value
pub const MAX_MOTIF_INDEX: usize = 4096;

/// Maximum confidence score for gene predictions
pub const MAX_CONFIDENCE_SCORE: f64 = 99.99;

/// Expected no-stop probability for genetic code analysis
pub const EXPECTED_NO_STOP_PROB: f64 = 0.953;

/// Test sequence repeat factor for creating longer sequences
pub const TEST_SEQUENCE_REPEAT_FACTOR: usize = 300;

/// Starting position for upstream analysis
pub const UPSTREAM_START_POS: usize = 1;

/// Ending position for upstream analysis
pub const UPSTREAM_END_POS: usize = 45;

/// GC content parameters for low AT frequency
pub const GC_LOW_AT_FREQ: f64 = 0.90;

/// GC content parameters for low GC frequency
pub const GC_LOW_GC_FREQ: f64 = 0.10;

/// GC content parameters for high AT frequency
pub const GC_HIGH_AT_FREQ: f64 = 0.10;

/// GC content parameters for high GC frequency
pub const GC_HIGH_GC_FREQ: f64 = 0.90;

// =============================================================================
// =============================================================================

/// Mask for lower bits in hamming distance calculations
pub const LOWER_BITS: u64 = 0x5555_5555_5555_5555;

/// Mask for upper bits in hamming distance calculations
pub const UPPER_BITS: u64 = 0xAAAA_AAAA_AAAA_AAAA;

/// Nucleotide lookup table for unpacking
pub const NUCLEOTIDE_LOOKUP: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Character lookup for sequence display
pub const NUCLEOTIDE_LETTERS: [char; 4] = ['A', 'G', 'C', 'T'];

// =============================================================================
// =============================================================================

/// Shine-Dalgarno motif descriptions for annotation
pub const RBS_DESCRIPTIONS: [(&str, &str); 28] = [
    ("None", "None"),
    ("GGA/GAG/AGG", "3-4bp"),
    ("3Base/5BMM", "13-15bp"),
    ("4Base/6BMM", "13-15bp"),
    ("AGxAG", "11-12bp"),
    ("AGxAG", "3-4bp"),
    ("GGA/GAG/AGG", "11-12bp"),
    ("GGxGG", "11-12bp"),
    ("GGxGG", "3-4bp"),
    ("AGxAG", "5-10bp"),
    ("AGGAG(G)/GGAGG", "13-15bp"),
    ("AGGA/GGAG/GAGG", "3-4bp"),
    ("AGGA/GGAG/GAGG", "11-12bp"),
    ("GGA/GAG/AGG", "5-10bp"),
    ("GGxGG", "5-10bp"),
    ("AGGA", "5-10bp"),
    ("GGAG/GAGG", "5-10bp"),
    ("AGxAGG/AGGxGG", "11-12bp"),
    ("AGxAGG/AGGxGG", "3-4bp"),
    ("AGxAGG/AGGxGG", "5-10bp"),
    ("AGGAG/GGAGG", "11-12bp"),
    ("AGGAG", "3-4bp"),
    ("AGGAG", "5-10bp"),
    ("GGAGG", "3-4bp"),
    ("GGAGG", "5-10bp"),
    ("AGGAGG", "11-12bp"),
    ("AGGAGG", "3-4bp"),
    ("AGGAGG", "5-10bp"),
];

/// Default start weight factor used in training and scoring
pub const DEFAULT_START_WEIGHT_FACTOR: f64 = 4.35;
