use crate::types::*;
use bio::io::fasta;
use std::fs::File;

/// Read sequence using rust-bio for FASTA files
/// Type alias to simplify the complex return type
pub type FastaRecord = (String, Option<String>, Vec<u8>);

pub fn read_fasta_sequences(filename: &str) -> Result<Vec<FastaRecord>, OrphosError> {
    let file = File::open(filename)?;
    let reader = fasta::Reader::new(file);
    let mut sequences = Vec::new();

    for result in reader.records() {
        let record = result.map_err(|e| OrphosError::ParseError(e.to_string()))?;
        let id = record.id().to_string();
        let description = record.desc().map(String::from);
        let seq = record.seq().to_vec();
        sequences.push((id, description, seq));
    }

    Ok(sequences)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_fasta_sequences_basic() {
        // Create a temporary FASTA file content
        let fasta_content = ">test_sequence\nATCG\nGCTA\n";

        // Write to temporary file for testing
        use std::env;
        use std::fs;
        let temp_dir = env::temp_dir();
        let temp_file = temp_dir.join("test_fasta.fa");
        fs::write(&temp_file, fasta_content).unwrap();

        let result = read_fasta_sequences(temp_file.to_str().unwrap());
        assert!(result.is_ok());

        let sequences = result.unwrap();
        assert_eq!(sequences.len(), 1);
        assert_eq!(sequences[0].0, "test_sequence");
        assert_eq!(sequences[0].2.len(), 8); // ATCGGCTA

        // Cleanup
        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_read_fasta_sequences_empty_file() {
        use std::env;
        use std::fs;
        let temp_dir = env::temp_dir();
        let temp_file = temp_dir.join("empty_fasta.fa");
        fs::write(&temp_file, "").unwrap();

        let result = read_fasta_sequences(temp_file.to_str().unwrap());
        assert!(result.is_ok());

        let sequences = result.unwrap();
        assert!(sequences.is_empty());

        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_read_fasta_sequences_multiple() {
        let fasta_content = ">seq1\nATCG\n>seq2\nGCTA\n>seq3\nTTAA\n";

        use std::env;
        use std::fs;
        let temp_dir = env::temp_dir();
        let temp_file = temp_dir.join("multi_fasta.fa");
        fs::write(&temp_file, fasta_content).unwrap();

        let result = read_fasta_sequences(temp_file.to_str().unwrap());
        assert!(result.is_ok());

        let sequences = result.unwrap();
        assert_eq!(sequences.len(), 3);
        assert_eq!(sequences[0].0, "seq1");
        assert_eq!(sequences[1].0, "seq2");
        assert_eq!(sequences[2].0, "seq3");

        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_read_fasta_sequences_with_description() {
        let fasta_content = ">seq1 This is a test sequence\nATCG\n>seq2\nGCTA\n";

        use std::env;
        use std::fs;
        let temp_dir = env::temp_dir();
        let temp_file = temp_dir.join("desc_fasta.fa");
        fs::write(&temp_file, fasta_content).unwrap();

        let result = read_fasta_sequences(temp_file.to_str().unwrap());
        assert!(result.is_ok());

        let sequences = result.unwrap();
        assert_eq!(sequences.len(), 2);
        assert_eq!(sequences[0].0, "seq1");
        assert_eq!(sequences[0].1, Some("This is a test sequence".to_string()));
        assert_eq!(sequences[1].1, None);

        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_read_fasta_sequences_file_not_found() {
        let result = read_fasta_sequences("nonexistent_file.fa");
        assert!(result.is_err());
        match result {
            Err(OrphosError::IoError(_)) => {}
            _ => panic!("Expected IoError for missing file"),
        }
    }

    #[test]
    fn test_read_fasta_sequences_invalid_format() {
        use std::env;
        use std::fs;
        let temp_dir = env::temp_dir();
        let temp_file = temp_dir.join("invalid_fasta.fa");
        // Create an invalid FASTA file (binary data)
        fs::write(&temp_file, vec![0x00, 0xFF, 0x80]).unwrap();

        let _result = read_fasta_sequences(temp_file.to_str().unwrap());
        let _ = fs::remove_file(temp_file);
    }

    #[test]
    fn test_fasta_record_type_alias() {
        // Test that the type alias works correctly
        let record: FastaRecord = (
            "test".to_string(),
            Some("desc".to_string()),
            vec![65, 84, 67, 71],
        );
        assert_eq!(record.0, "test");
        assert_eq!(record.1, Some("desc".to_string()));
        assert_eq!(record.2, vec![65, 84, 67, 71]);
    }
}
