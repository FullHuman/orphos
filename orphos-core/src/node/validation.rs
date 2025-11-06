use crate::{
    constants::MINIMUM_EDGE_GENE_LENGTH,
    sequence::{is_atg, is_gtg, is_start, is_ttg},
    types::{CodonType, Mask, Training},
};

/// Check if position contains a start codon and return its type
#[inline]
pub fn check_start_codon(
    sequence: &[u8],
    position: usize,
    training: &Training,
) -> Option<CodonType> {
    if !is_start(sequence, position, training) {
        return None;
    }

    if is_atg(sequence, position) {
        Some(CodonType::Atg)
    } else if is_gtg(sequence, position) {
        Some(CodonType::Gtg)
    } else if is_ttg(sequence, position) {
        Some(CodonType::Ttg)
    } else {
        None
    }
}

/// Check if gene is valid (length and mask constraints)
pub fn is_valid_gene(
    start_position: usize,
    stop_position: usize,
    minimum_gene_distance: usize,
    masks: &[Mask],
) -> bool {
    (stop_position + 3).saturating_sub(start_position) >= minimum_gene_distance
        && !crosses_mask(start_position, stop_position, masks)
}

/// Check if gene is valid for reverse strand
pub fn is_valid_reverse_gene(
    start_position: usize,
    stop_position: usize,
    minimum_gene_distance: usize,
    sequence_length: usize,
    masks: &[Mask],
) -> bool {
    (stop_position + 3).saturating_sub(start_position) >= minimum_gene_distance
        && !crosses_mask(
            sequence_length - stop_position - 1,
            sequence_length - start_position - 1,
            masks,
        )
}

/// Check if position qualifies as edge gene
pub fn is_edge_gene(position: usize, stop: usize, closed: bool, masks: &[Mask]) -> bool {
    position <= 2
        && !closed
        && stop.saturating_sub(position) > MINIMUM_EDGE_GENE_LENGTH
        && !crosses_mask(position, stop, masks)
}

/// Check if position qualifies as reverse edge gene
pub fn is_reverse_edge_gene(
    position: usize,
    stop: usize,
    sequence_length: usize,
    closed: bool,
    masks: &[Mask],
) -> bool {
    position <= 2
        && !closed
        && stop.saturating_sub(position) > MINIMUM_EDGE_GENE_LENGTH
        && !crosses_mask(
            sequence_length - stop - 1,
            sequence_length - position - 1,
            masks,
        )
}

/// Check if a gene boundary crosses any mask
fn crosses_mask(start: usize, end: usize, masks: &[Mask]) -> bool {
    masks
        .iter()
        .any(|mask| !(end < mask.begin || start > mask.end))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::encode_sequence;
    use crate::types::{Mask, Training};

    fn get_encoded_sequence(input: &[u8]) -> Vec<u8> {
        let sequence_length = input.len();
        let mut seq = vec![0u8; (sequence_length * 2).div_ceil(8)]; // 2 bits per nucleotide
        let mut unknown_sequence = vec![0u8; sequence_length.div_ceil(8)]; // 1 bit per nucleotide for unknowns
        let mut masks = Vec::new();
        let _ = encode_sequence(input, &mut seq, &mut unknown_sequence, &mut masks, false).unwrap();
        seq
    }

    #[test]
    fn test_check_start_codon_atg() {
        let sequence = get_encoded_sequence(b"ATGCCC");
        let training = Training::default();
        assert_eq!(
            check_start_codon(&sequence, 0, &training),
            Some(CodonType::Atg)
        );
    }

    #[test]
    fn test_check_start_codon_gtg() {
        let sequence = get_encoded_sequence(b"GTGCCC");
        let training = Training::default();
        assert_eq!(
            check_start_codon(&sequence, 0, &training),
            Some(CodonType::Gtg)
        );
    }

    #[test]
    fn test_check_start_codon_ttg() {
        let sequence = get_encoded_sequence(b"TTGCCC");
        let training = Training::default();
        assert_eq!(
            check_start_codon(&sequence, 0, &training),
            Some(CodonType::Ttg)
        );
    }

    #[test]
    fn test_check_start_codon_none() {
        let sequence = get_encoded_sequence(b"CCCGGG");
        let training = Training::default();
        assert_eq!(check_start_codon(&sequence, 0, &training), None);
    }

    #[test]
    fn test_is_valid_gene_valid() {
        let masks = vec![];
        assert!(is_valid_gene(0, 100, 90, &masks));
    }

    #[test]
    fn test_is_valid_gene_too_short() {
        let masks = vec![];
        assert!(!is_valid_gene(0, 50, 90, &masks));
    }

    #[test]
    fn test_is_valid_gene_crosses_mask() {
        let masks = vec![Mask { begin: 10, end: 20 }];
        assert!(!is_valid_gene(5, 100, 90, &masks));
    }

    #[test]
    fn test_is_valid_reverse_gene_valid() {
        let masks = vec![];
        assert!(is_valid_reverse_gene(0, 100, 90, 1000, &masks));
    }

    #[test]
    fn test_is_valid_reverse_gene_too_short() {
        let masks = vec![];
        assert!(!is_valid_reverse_gene(0, 50, 90, 1000, &masks));
    }

    #[test]
    fn test_is_edge_gene_valid() {
        let masks = vec![];
        assert!(is_edge_gene(1, 200, false, &masks));
    }

    #[test]
    fn test_is_edge_gene_position_too_high() {
        let masks = vec![];
        assert!(!is_edge_gene(5, 200, false, &masks));
    }

    #[test]
    fn test_is_edge_gene_closed() {
        let masks = vec![];
        assert!(!is_edge_gene(1, 200, true, &masks));
    }

    #[test]
    fn test_is_edge_gene_too_short() {
        let masks = vec![];
        assert!(!is_edge_gene(1, 50, false, &masks));
    }

    #[test]
    fn test_is_reverse_edge_gene_valid() {
        let masks = vec![];
        assert!(is_reverse_edge_gene(1, 200, 1000, false, &masks));
    }

    #[test]
    fn test_is_reverse_edge_gene_closed() {
        let masks = vec![];
        assert!(!is_reverse_edge_gene(1, 200, 1000, true, &masks));
    }

    #[test]
    fn test_crosses_mask_no_overlap() {
        let masks = vec![Mask {
            begin: 50,
            end: 100,
        }];
        assert!(!crosses_mask(10, 30, &masks));
        assert!(!crosses_mask(120, 150, &masks));
    }

    #[test]
    fn test_crosses_mask_overlap() {
        let masks = vec![Mask {
            begin: 50,
            end: 100,
        }];
        assert!(crosses_mask(40, 60, &masks));
        assert!(crosses_mask(80, 120, &masks));
        assert!(crosses_mask(30, 150, &masks));
    }

    #[test]
    fn test_crosses_mask_multiple_masks() {
        let masks = vec![Mask { begin: 10, end: 20 }, Mask { begin: 50, end: 60 }];
        assert!(crosses_mask(5, 15, &masks));
        assert!(crosses_mask(55, 65, &masks));
        assert!(!crosses_mask(25, 45, &masks));
    }
}
