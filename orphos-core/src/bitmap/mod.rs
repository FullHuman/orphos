#[inline]
const fn calculate_bit_position(bit_index: usize) -> (usize, u8) {
    (bit_index >> 3, 1 << (bit_index & 0x07))
}

/// Test if a bit is set at the given index
pub fn test_bit(bitmap: &[u8], bit_index: usize) -> bool {
    let (byte_index, bit_mask) = calculate_bit_position(bit_index);
    (bitmap[byte_index] & bit_mask) != 0
}

/// Set a bit to 1 at the given index
pub fn set_bit(bitmap: &mut [u8], bit_index: usize) {
    let (byte_index, bit_mask) = calculate_bit_position(bit_index);
    bitmap[byte_index] |= bit_mask;
}

/// Clear a bit (set it to 0) at the given index
pub fn clear_bit(bitmap: &mut [u8], bit_index: usize) {
    let (byte_index, bit_mask) = calculate_bit_position(bit_index);
    bitmap[byte_index] &= !bit_mask;
}

/// Flip a bit's value 0->1 or 1->0 at the given index
pub fn toggle_bit(bitmap: &mut [u8], bit_index: usize) {
    let (byte_index, bit_mask) = calculate_bit_position(bit_index);
    bitmap[byte_index] ^= bit_mask;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_bit_position() {
        assert_eq!(calculate_bit_position(0), (0, 1));
        assert_eq!(calculate_bit_position(1), (0, 2));
        assert_eq!(calculate_bit_position(7), (0, 128));
        assert_eq!(calculate_bit_position(8), (1, 1));
        assert_eq!(calculate_bit_position(15), (1, 128));
        assert_eq!(calculate_bit_position(16), (2, 1));
    }

    #[test]
    fn test_test_bit() {
        let bitmap = [0b10101010, 0b01010101];

        // Test bits in first byte
        assert!(!test_bit(&bitmap, 0));
        assert!(test_bit(&bitmap, 1));
        assert!(!test_bit(&bitmap, 2));
        assert!(test_bit(&bitmap, 3));
        assert!(!test_bit(&bitmap, 4));
        assert!(test_bit(&bitmap, 5));
        assert!(!test_bit(&bitmap, 6));
        assert!(test_bit(&bitmap, 7));

        // Test bits in second byte
        assert!(test_bit(&bitmap, 8));
        assert!(!test_bit(&bitmap, 9));
        assert!(test_bit(&bitmap, 10));
        assert!(!test_bit(&bitmap, 11));
    }

    #[test]
    fn test_set_bit() {
        let mut bitmap = [0u8; 2];

        set_bit(&mut bitmap, 0);
        assert_eq!(bitmap[0], 0b00000001);

        set_bit(&mut bitmap, 3);
        assert_eq!(bitmap[0], 0b00001001);

        set_bit(&mut bitmap, 8);
        assert_eq!(bitmap[1], 0b00000001);

        set_bit(&mut bitmap, 0);
        assert_eq!(bitmap[0], 0b00001001);
    }

    #[test]
    fn test_clear_bit() {
        let mut bitmap = [0b11111111, 0b11111111];

        clear_bit(&mut bitmap, 0);
        assert_eq!(bitmap[0], 0b11111110);

        clear_bit(&mut bitmap, 3);
        assert_eq!(bitmap[0], 0b11110110);

        clear_bit(&mut bitmap, 8);
        assert_eq!(bitmap[1], 0b11111110);

        clear_bit(&mut bitmap, 0);
        assert_eq!(bitmap[0], 0b11110110);
    }

    #[test]
    fn test_toggle_bit() {
        let mut bitmap = [0b10101010, 0b01010101];

        // Toggle bit from 0 to 1
        toggle_bit(&mut bitmap, 0);
        assert_eq!(bitmap[0], 0b10101011);

        // Toggle bit from 1 to 0
        toggle_bit(&mut bitmap, 1);
        assert_eq!(bitmap[0], 0b10101001);

        // Toggle in second byte
        toggle_bit(&mut bitmap, 8);
        assert_eq!(bitmap[1], 0b01010100);

        toggle_bit(&mut bitmap, 8);
        assert_eq!(bitmap[1], 0b01010101);
    }

    #[test]
    fn test_bit_operations_consistency() {
        let mut bitmap = [0u8; 4];

        for i in 0..32 {
            assert!(!test_bit(&bitmap, i));

            // Set bit and verify
            set_bit(&mut bitmap, i);
            assert!(test_bit(&bitmap, i));

            toggle_bit(&mut bitmap, i);
            assert!(!test_bit(&bitmap, i));

            toggle_bit(&mut bitmap, i);
            assert!(test_bit(&bitmap, i));

            clear_bit(&mut bitmap, i);
            assert!(!test_bit(&bitmap, i));
        }
    }

    #[test]
    fn test_edge_cases() {
        let mut bitmap = [0u8; 1];

        // Test all 8 bits in a single byte
        for i in 0..8 {
            set_bit(&mut bitmap, i);
        }
        assert_eq!(bitmap[0], 0b11111111);

        for i in 0..8 {
            clear_bit(&mut bitmap, i);
        }
        assert_eq!(bitmap[0], 0b00000000);
    }
}
