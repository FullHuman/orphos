use crate::{
    sequence::{
        create_reverse_complement_sequence, encode_sequence, encode_sequence_simd_wide_packed,
    },
    types::Mask,
};

#[derive(Debug)]
pub struct EncodedSequence {
    pub forward_sequence: Vec<u8>,
    pub reverse_complement_sequence: Vec<u8>,
    pub unknown_sequence: Vec<u8>,
    pub masks: Vec<Mask>,
    pub gc_content: f64,
    pub sequence_length: usize,
}

impl EncodedSequence {
    pub fn with_masking(sequence: &[u8]) -> Self {
        let nucleotide_length = sequence.len();
        let mut forward_sequence: Vec<u8> = vec![0; (nucleotide_length * 2).div_ceil(8)];
        let mut unknown_sequence: Vec<u8> = vec![0; nucleotide_length.div_ceil(8)];

        let mut masks = vec![];

        let gc_content = encode_sequence(
            sequence,
            &mut forward_sequence,
            &mut unknown_sequence,
            &mut masks,
            true,
        )
        .unwrap();
        let reverse_complement_sequence = create_reverse_complement_sequence(
            &forward_sequence,
            &unknown_sequence,
            nucleotide_length,
        );

        Self {
            forward_sequence,
            reverse_complement_sequence,
            unknown_sequence,
            masks,
            gc_content,
            sequence_length: nucleotide_length,
        }
    }

    pub fn without_masking(sequence: &[u8]) -> Self {
        let nucleotide_length = sequence.len();
        let mut forward_sequence: Vec<u8> = vec![0; (nucleotide_length * 2).div_ceil(8)];
        let mut unknown_sequence: Vec<u8> = vec![0; nucleotide_length.div_ceil(8)];
        let masks = vec![];

        // let gc_content = encode_sequence(sequence, &mut forward_sequence, &mut unknown_sequence, &mut masks, false).unwrap();
        let gc_content = encode_sequence_simd_wide_packed(
            sequence,
            &mut forward_sequence,
            &mut unknown_sequence,
        )
        .unwrap();
        let reverse_complement_sequence = create_reverse_complement_sequence(
            &forward_sequence,
            &unknown_sequence,
            nucleotide_length,
        );

        Self {
            forward_sequence,
            reverse_complement_sequence,
            unknown_sequence,
            masks,
            gc_content,
            sequence_length: nucleotide_length,
        }
    }
}
