use std::io::Write;

use bio::bio_types::strand::Strand;

use crate::{OrphosError, constants::VERSION, results::OrphosResults};

/// Write results in GFF format
pub fn write_gff_format<W: Write>(
    writer: &mut W,
    results: &OrphosResults,
) -> Result<(), OrphosError> {
    writeln!(writer, "##gff-version  3")?;
    // Add header comments like original Orphos
    if let Some(desc) = &results.sequence_info.description {
        writeln!(
            writer,
            "# Sequence Data: seqnum=1;seqlen={};seqhdr=\"{} {}\"",
            results.sequence_info.length, results.sequence_info.header, desc
        )?;
    } else {
        writeln!(
            writer,
            "# Sequence Data: seqnum=1;seqlen={};seqhdr=\"{}\"",
            results.sequence_info.length, results.sequence_info.header
        )?;
    }
    writeln!(
        writer,
        "# Model Data: version=Orphos.v{};run_type={};model=\"{}\";gc_cont={:.2};transl_table={};uses_sd={}",
        VERSION,
        if results.metagenomic_model.is_some() {
            "Meta"
        } else {
            "Single"
        },
        results.metagenomic_model.as_deref().unwrap_or("Ab initio"),
        results.sequence_info.gc_content * 100.0,
        results.training_used.translation_table,
        if results.training_used.uses_shine_dalgarno {
            1
        } else {
            0
        }
    )?;

    for gene in results.genes.iter() {
        let strand_char = match gene.coordinates.strand {
            Strand::Forward => '+',
            Strand::Reverse => '-',
            Strand::Unknown => '.',
        };

        writeln!(
            writer,
            "{}\tOrphos_v{}\tCDS\t{}\t{}\t{:.1}\t{}\t0\t{};{};",
            results.sequence_info.header,
            VERSION,
            gene.coordinates.begin,
            gene.coordinates.end,
            gene.score.coding_score + gene.score.start_score,
            strand_char,
            gene.annotation,
            gene.score
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::results::{OrphosResults, SequenceInfo};
    use crate::types::*;
    use bio::bio_types::strand::Strand;

    fn create_test_results() -> OrphosResults {
        let sequence_info = SequenceInfo {
            header: "test_sequence".to_string(),
            description: Some("Test sequence for unit tests".to_string()),
            length: 1000,
            gc_content: 0.45,
            num_genes: 2,
        };

        let training = Training {
            translation_table: 11,
            uses_shine_dalgarno: true,
            ..Default::default()
        };

        let gene1 = Gene {
            coordinates: GeneCoordinates {
                begin: 100,
                end: 300,
                strand: Strand::Forward,
                start_index: 100,
                stop_index: 300,
            },
            score: GeneScore {
                confidence: 95.5,
                total_score: 12.5,
                coding_score: 8.2,
                start_score: 4.3,
                ribosome_binding_score: 3.1,
                upstream_score: 2.8,
                type_score: 1.9,
            },
            annotation: GeneAnnotation::new(
                "gene_1".to_string(),
                false,
                false,
                StartType::Atg,
                0.42,
            )
            .with_rbs("AGGAG".to_string(), "5".to_string()),
        };

        let gene2 = Gene {
            coordinates: GeneCoordinates {
                begin: 400,
                end: 600,
                strand: Strand::Reverse,
                start_index: 600,
                stop_index: 400,
            },
            score: GeneScore {
                confidence: 87.3,
                total_score: 10.8,
                coding_score: 7.1,
                start_score: 3.7,
                ribosome_binding_score: 2.9,
                upstream_score: 2.1,
                type_score: 1.2,
            },
            annotation: GeneAnnotation::new(
                "gene_2".to_string(),
                true,
                false,
                StartType::Gtg,
                0.38,
            ),
        };

        OrphosResults {
            sequence_info,
            training_used: training,
            genes: vec![gene1, gene2],
            metagenomic_model: None,
        }
    }

    #[test]
    fn test_gff_header_with_description() {
        let results = create_test_results();
        let mut output = Vec::new();

        write_gff_format(&mut output, &results).unwrap();
        let gff_output = String::from_utf8(output).unwrap();

        assert!(gff_output.starts_with("##gff-version  3\n"));
        assert!(gff_output.contains("# Sequence Data: seqnum=1;seqlen=1000;seqhdr=\"test_sequence Test sequence for unit tests\""));
        assert!(gff_output.contains("# Model Data: version=Orphos.v"));
        assert!(gff_output.contains("run_type=Single"));
        assert!(gff_output.contains("model=\"Ab initio\""));
        assert!(gff_output.contains("gc_cont=45.00"));
        assert!(gff_output.contains("transl_table=11"));
        assert!(gff_output.contains("uses_sd=1"));
    }

    #[test]
    fn test_gff_header_without_description() {
        let mut results = create_test_results();
        results.sequence_info.description = None;
        let mut output = Vec::new();

        write_gff_format(&mut output, &results).unwrap();
        let gff_output = String::from_utf8(output).unwrap();

        assert!(
            gff_output.contains("# Sequence Data: seqnum=1;seqlen=1000;seqhdr=\"test_sequence\"")
        );
        assert!(!gff_output.contains("Test sequence for unit tests"));
    }

    #[test]
    fn test_gff_metagenomic_model() {
        let mut results = create_test_results();
        results.metagenomic_model = Some("Meta_model_v1".to_string());
        let mut output = Vec::new();

        write_gff_format(&mut output, &results).unwrap();
        let gff_output = String::from_utf8(output).unwrap();

        assert!(gff_output.contains("run_type=Meta"));
        assert!(gff_output.contains("model=\"Meta_model_v1\""));
    }

    #[test]
    fn test_gff_gene_entries() {
        let results = create_test_results();
        let mut output = Vec::new();

        write_gff_format(&mut output, &results).unwrap();
        let gff_output = String::from_utf8(output).unwrap();

        // Check forward strand gene
        assert!(gff_output.contains("test_sequence\tOrphos_v"));
        assert!(gff_output.contains("\tCDS\t100\t300\t12.5\t+\t0\t"));

        // Check reverse strand gene
        assert!(gff_output.contains("\tCDS\t400\t600\t10.8\t-\t0\t"));

        // Check that gene annotations are included
        assert!(gff_output.contains("ID=gene_1"));
        assert!(gff_output.contains("ID=gene_2"));
    }

    #[test]
    fn test_gff_strand_characters() {
        let mut results = create_test_results();
        results.genes[0].coordinates.strand = Strand::Forward;
        results.genes[1].coordinates.strand = Strand::Reverse;

        // Add unknown strand gene
        let mut unknown_gene = results.genes[0].clone();
        unknown_gene.coordinates.strand = Strand::Unknown;
        unknown_gene.annotation.identifier = "gene_3".to_string();
        results.genes.push(unknown_gene);

        let mut output = Vec::new();
        write_gff_format(&mut output, &results).unwrap();
        let gff_output = String::from_utf8(output).unwrap();

        assert!(gff_output.contains("\t+\t"));
        assert!(gff_output.contains("\t-\t"));
        assert!(gff_output.contains("\t.\t"));
    }

    #[test]
    fn test_gff_no_shine_dalgarno() {
        let mut results = create_test_results();
        results.training_used.uses_shine_dalgarno = false;
        let mut output = Vec::new();

        write_gff_format(&mut output, &results).unwrap();
        let gff_output = String::from_utf8(output).unwrap();

        assert!(gff_output.contains("uses_sd=0"));
    }

    #[test]
    fn test_gff_empty_genes() {
        let mut results = create_test_results();
        results.genes.clear();
        let mut output = Vec::new();

        write_gff_format(&mut output, &results).unwrap();
        let gff_output = String::from_utf8(output).unwrap();

        // Should still have headers but no gene entries
        assert!(gff_output.starts_with("##gff-version  3\n"));
        assert!(gff_output.contains("# Sequence Data:"));
        assert!(gff_output.contains("# Model Data:"));
        assert!(!gff_output.contains("\tCDS\t"));
    }
}
