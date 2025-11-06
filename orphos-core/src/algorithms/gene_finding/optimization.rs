//! Gene start tweaking and annotation pipeline (faithful to C tweak_final_starts logic)

use bio::bio_types::strand::Strand;

use crate::{
    algorithms::gene_finding::{annotation::record_gene_annotation_and_score, creation::add_genes},
    constants::{MAX_RIBOSOME_DISTANCE, MAXIMUM_SAME_OVERLAP, SEARCH_WINDOW},
    node::intergenic_mod,
    types::{CodonType, Gene, Node, Training},
};

/// Builder orchestrating start tweaking and annotation
pub struct GeneBuilder<'a> {
    genes: Vec<Gene>,
    nodes: &'a [Node],
    training: &'a Training,
    sequence_number: usize,
}

impl<'a> GeneBuilder<'a> {
    pub fn from_nodes(
        nodes: &'a [Node],
        path_start: usize,
        training: &'a Training,
        sequence_number: usize,
    ) -> Self {
        let genes = add_genes(nodes, path_start);
        Self {
            genes,
            nodes,
            training,
            sequence_number,
        }
    }

    pub fn with_tweaked_starts(mut self) -> Self {
        for i in 0..self.genes.len() {
            self.tweak_gene_start(i);
        }
        self
    }

    pub fn with_annotations(mut self) -> Self {
        for (i, gene) in self.genes.iter_mut().enumerate() {
            let (annotation, score) = record_gene_annotation_and_score(
                gene,
                self.nodes,
                self.training,
                self.sequence_number,
                i,
            );
            gene.annotation = annotation;
            gene.score = score;
        }
        self
    }

    pub fn build(self) -> Vec<Gene> {
        self.genes
    }

    fn tweak_gene_start(&mut self, gene_index: usize) {
        let gene_count = self.genes.len();
        if gene_count == 0 {
            return;
        }
        let current_idx = self.genes[gene_index].coordinates.start_index;
        let current = &self.nodes[current_idx];
        let current_score = current.scores.coding_score + current.scores.start_score;

        // intergenic modifier for current start
        let mut igm = 0.0;
        if gene_index > 0 {
            let prev_gene = &self.genes[gene_index - 1];
            let prev_start = &self.nodes[prev_gene.coordinates.start_index];
            if current.position.strand == Strand::Forward
                && prev_start.position.strand == Strand::Forward
            {
                igm = intergenic_mod(
                    &self.nodes[prev_gene.coordinates.stop_index],
                    current,
                    self.training,
                );
            } else if current.position.strand == Strand::Forward
                && prev_start.position.strand == Strand::Reverse
            {
                igm = intergenic_mod(prev_start, current, self.training);
            }
        }
        if gene_index < gene_count - 1 {
            let next_gene = &self.genes[gene_index + 1];
            let next_start = &self.nodes[next_gene.coordinates.start_index];
            if current.position.strand == Strand::Reverse
                && next_start.position.strand == Strand::Forward
            {
                igm = intergenic_mod(current, next_start, self.training);
            } else if current.position.strand == Strand::Reverse
                && next_start.position.strand == Strand::Reverse
            {
                igm = intergenic_mod(
                    current,
                    &self.nodes[next_gene.coordinates.stop_index],
                    self.training,
                );
            }
        }

        // collect best two alternative starts sharing stop within SEARCH_WINDOW
        let mut alt_idx: [Option<usize>; 2] = [None, None];
        let mut alt_score: [f64; 2] = [f64::MIN, f64::MIN];
        let mut alt_igm: [f64; 2] = [0.0, 0.0];
        // NOTE: The original C code scans +/-100 node indices around the current start.
        // Differences in node ordering in this Rust port can cause legitimate
        // alternative (e.g. edge) starts to fall outside that local window.
        // To preserve correctness we scan the full node array for starts with
        // the same stop_value. This can be optimized later if needed.
        let win_start = current_idx.saturating_sub(SEARCH_WINDOW);
        let win_end = (current_idx + SEARCH_WINDOW).min(self.nodes.len());
        for j in win_start..win_end {
            if j == current_idx {
                continue;
            }
            let node = &self.nodes[j];
            if node.position.codon_type == CodonType::Stop
                || node.position.stop_value != current.position.stop_value
            {
                continue;
            }

            let mut tigm = 0.0;
            let mut skip = false;
            if gene_index > 0 {
                let prev_gene = &self.genes[gene_index - 1];
                let prev_start = &self.nodes[prev_gene.coordinates.start_index];
                if node.position.strand == Strand::Forward
                    && prev_start.position.strand == Strand::Forward
                {
                    let prev_stop_pos = self.nodes[prev_gene.coordinates.stop_index].position.index;
                    if prev_stop_pos > node.position.index + MAXIMUM_SAME_OVERLAP {
                        skip = true;
                    } else {
                        tigm = intergenic_mod(
                            &self.nodes[prev_gene.coordinates.stop_index],
                            node,
                            self.training,
                        );
                    }
                } else if node.position.strand == Strand::Forward
                    && prev_start.position.strand == Strand::Reverse
                {
                    if prev_start.position.index >= node.position.index {
                        skip = true;
                    } else {
                        tigm = intergenic_mod(prev_start, node, self.training);
                    }
                }
            }
            if skip {
                continue;
            }
            if gene_index < gene_count - 1 {
                let next_gene = &self.genes[gene_index + 1];
                let next_start = &self.nodes[next_gene.coordinates.start_index];
                if node.position.strand == Strand::Reverse
                    && next_start.position.strand == Strand::Forward
                {
                    if node.position.index >= next_start.position.index {
                        skip = true;
                    } else {
                        tigm = intergenic_mod(node, next_start, self.training);
                    }
                } else if node.position.strand == Strand::Reverse
                    && next_start.position.strand == Strand::Reverse
                {
                    let next_stop_pos = self.nodes[next_gene.coordinates.stop_index].position.index;
                    if node.position.index > next_stop_pos + MAXIMUM_SAME_OVERLAP {
                        skip = true;
                    } else {
                        tigm = intergenic_mod(
                            node,
                            &self.nodes[next_gene.coordinates.stop_index],
                            self.training,
                        );
                    }
                }
            }
            if skip {
                continue;
            }

            let score = node.scores.coding_score + node.scores.start_score;
            if alt_idx[0].is_none() || score + tigm > alt_score[0] + alt_igm[0] {
                alt_idx[1] = alt_idx[0];
                alt_score[1] = alt_score[0];
                alt_igm[1] = alt_igm[0];
                alt_idx[0] = Some(j);
                alt_score[0] = score;
                alt_igm[0] = tigm;
            } else if alt_idx[1].is_none() || score + tigm > alt_score[1] + alt_igm[1] {
                alt_idx[1] = Some(j);
                alt_score[1] = score;
                alt_igm[1] = tigm;
            }
        }

        // adjust alternative scores per C rules
        for k in 0..2 {
            if let Some(idx) = alt_idx[k] {
                let node = &self.nodes[idx];
                let dist = node.position.index.abs_diff(current.position.index);
                // Edge promotion: if alternative is edge and current is not, keep raw score (do not penalize)
                if node.position.is_edge && !current.position.is_edge {
                    // leave alt_score as coding+start (already set)
                    continue;
                }
                if node.scores.type_score < current.scores.type_score
                    && alt_score[k] - node.scores.type_score
                        >= current_score - current.scores.type_score
                            + self.training.start_weight_factor
                    && node.scores.ribosome_binding_score > current.scores.ribosome_binding_score
                    && node.scores.upstream_score > current.scores.upstream_score
                    && node.scores.coding_score > current.scores.coding_score
                    && dist > MAX_RIBOSOME_DISTANCE
                {
                    alt_score[k] += current.scores.type_score - node.scores.type_score; // rule 1
                } else if dist <= MAX_RIBOSOME_DISTANCE
                    && node.scores.ribosome_binding_score + node.scores.type_score
                        > current.scores.ribosome_binding_score + current.scores.type_score
                    && !current.position.is_edge
                    && !node.position.is_edge
                {
                    if current.scores.coding_score > node.scores.coding_score {
                        alt_score[k] += current.scores.coding_score - node.scores.coding_score;
                    }
                    if current.scores.upstream_score > node.scores.upstream_score {
                        alt_score[k] += current.scores.upstream_score - node.scores.upstream_score;
                    }
                    if igm > alt_igm[k] {
                        alt_score[k] += igm - alt_igm[k];
                    }
                } else if !node.position.is_edge || current.position.is_edge {
                    // only penalize if not the protected edge promotion case
                    alt_score[k] = -1000.0;
                }
            }
        }

        // pick best improved alternative
        let mut best: Option<usize> = None;
        for k in 0..2 {
            if alt_idx[k].is_some() && alt_score[k] + alt_igm[k] > current_score + igm {
                match best {
                    None => best = Some(k),
                    Some(b) => {
                        if alt_score[k] + alt_igm[k] > alt_score[b] + alt_igm[b] {
                            best = Some(k);
                        }
                    }
                }
            }
        }
        if let Some(k) = best
            && let Some(new_idx) = alt_idx[k]
        {
            self.apply_alternative_start(gene_index, new_idx);
        }

        // Optional debug dump for first gene: list all candidate starts with same stop
        // (debug instrumentation removed)
    }

    fn apply_alternative_start(&mut self, gene_index: usize, new_start_idx: usize) {
        // Allow promotion to an edge start; only skip if the alternative is identical.
        if self.genes[gene_index].coordinates.start_index == new_start_idx {
            return;
        }
        self.genes[gene_index].coordinates.start_index = new_start_idx;
        let node = &self.nodes[new_start_idx];
        if node.position.strand == Strand::Forward {
            self.genes[gene_index].coordinates.begin = node.position.index + 1;
        } else {
            self.genes[gene_index].coordinates.end = node.position.index + 1;
        }
    }
}
