#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MetagenomicBin {
    /// Zero-based bin identifier (0..49)
    pub id: usize,
    /// Canonical organism label used in C implementation
    pub name: &'static str,
    /// Domain: "A" = Archaea, "B" = Bacteria (from C metadata)
    pub domain: &'static str,
    /// Representative genome GC percentage (0-100 scale)
    pub gc_percent: f64,
}

/// The 50 metagenomic bins (descriptors only).
pub const fn bins() -> &'static [MetagenomicBin] {
    &BINS
}

const BINS: [MetagenomicBin; 50] = [
    MetagenomicBin {
        id: 0,
        name: "Mycoplasma_bovis_PG45",
        domain: "B",
        gc_percent: 29.31,
    },
    MetagenomicBin {
        id: 1,
        name: "Mycoplasma_pneumoniae_M129",
        domain: "B",
        gc_percent: 40.01,
    },
    MetagenomicBin {
        id: 2,
        name: "Mycoplasma_suis_Illinois",
        domain: "B",
        gc_percent: 31.08,
    },
    MetagenomicBin {
        id: 3,
        name: "Aeropyrum_pernix_K1",
        domain: "A",
        gc_percent: 56.31,
    },
    MetagenomicBin {
        id: 4,
        name: "Akkermansia_muciniphila_ATCC_BAA_835",
        domain: "B",
        gc_percent: 55.76,
    },
    MetagenomicBin {
        id: 5,
        name: "Anaplasma_marginale_Maries",
        domain: "B",
        gc_percent: 49.76,
    },
    MetagenomicBin {
        id: 6,
        name: "Anaplasma_phagocytophilum_HZ",
        domain: "B",
        gc_percent: 41.64,
    },
    MetagenomicBin {
        id: 7,
        name: "Archaeoglobus_fulgidus_DSM_4304",
        domain: "A",
        gc_percent: 48.58,
    },
    MetagenomicBin {
        id: 8,
        name: "Bacteroides_fragilis_NCTC_9343",
        domain: "B",
        gc_percent: 43.19,
    },
    MetagenomicBin {
        id: 9,
        name: "Brucella_canis_ATCC_23365",
        domain: "B",
        gc_percent: 57.21,
    },
    MetagenomicBin {
        id: 10,
        name: "Burkholderia_rhizoxinica_HKI_454",
        domain: "B",
        gc_percent: 59.70,
    },
    MetagenomicBin {
        id: 11,
        name: "Candidatus_Amoebophilus_asiaticus_5a2",
        domain: "B",
        gc_percent: 35.05,
    },
    MetagenomicBin {
        id: 12,
        name: "Candidatus_Korarchaeum_cryptofilum_OPF8",
        domain: "A",
        gc_percent: 49.00,
    },
    MetagenomicBin {
        id: 13,
        name: "Catenulispora_acidiphila_DSM_44928",
        domain: "B",
        gc_percent: 69.77,
    },
    MetagenomicBin {
        id: 14,
        name: "Cenarchaeum_symbiosum_B",
        domain: "A",
        gc_percent: 57.19,
    },
    MetagenomicBin {
        id: 15,
        name: "Chlorobium_phaeobacteroides_BS1",
        domain: "B",
        gc_percent: 48.93,
    },
    MetagenomicBin {
        id: 16,
        name: "Chlorobium_tepidum_TLS",
        domain: "B",
        gc_percent: 56.53,
    },
    MetagenomicBin {
        id: 17,
        name: "Desulfotomaculum_acetoxidans_DSM_771",
        domain: "B",
        gc_percent: 41.55,
    },
    MetagenomicBin {
        id: 18,
        name: "Desulfurococcus_kamchatkensis_1221n",
        domain: "B",
        gc_percent: 45.34,
    },
    MetagenomicBin {
        id: 19,
        name: "Erythrobacter_litoralis_HTCC2594",
        domain: "B",
        gc_percent: 63.07,
    },
    MetagenomicBin {
        id: 20,
        name: "Escherichia_coli_UMN026",
        domain: "B",
        gc_percent: 50.72,
    },
    MetagenomicBin {
        id: 21,
        name: "Haloquadratum_walsbyi_DSM_16790",
        domain: "A",
        gc_percent: 47.86,
    },
    MetagenomicBin {
        id: 22,
        name: "Halorubrum_lacusprofundi_ATCC_49239",
        domain: "A",
        gc_percent: 57.14,
    },
    MetagenomicBin {
        id: 23,
        name: "Hyperthermus_butylicus_DSM_5456",
        domain: "A",
        gc_percent: 53.74,
    },
    MetagenomicBin {
        id: 24,
        name: "Ignisphaera_aggregans_DSM_17230",
        domain: "A",
        gc_percent: 35.69,
    },
    MetagenomicBin {
        id: 25,
        name: "Marinobacter_aquaeolei_VT8",
        domain: "B",
        gc_percent: 57.27,
    },
    MetagenomicBin {
        id: 26,
        name: "Methanopyrus_kandleri_AV19",
        domain: "A",
        gc_percent: 61.16,
    },
    MetagenomicBin {
        id: 27,
        name: "Methanosphaerula_palustris_E1_9c",
        domain: "A",
        gc_percent: 55.35,
    },
    MetagenomicBin {
        id: 28,
        name: "Methanothermobacter_thermautotrophicus_Delta_H",
        domain: "B",
        gc_percent: 49.54,
    },
    MetagenomicBin {
        id: 29,
        name: "Methylacidiphilum_infernorum_V4",
        domain: "B",
        gc_percent: 45.48,
    },
    MetagenomicBin {
        id: 30,
        name: "Mycobacterium_leprae_TN",
        domain: "B",
        gc_percent: 57.80,
    },
    MetagenomicBin {
        id: 31,
        name: "Natrialba_magadii_ATCC_43099",
        domain: "A",
        gc_percent: 61.42,
    },
    MetagenomicBin {
        id: 32,
        name: "Orientia_tsutsugamushi_Boryong",
        domain: "B",
        gc_percent: 30.53,
    },
    MetagenomicBin {
        id: 33,
        name: "Pelotomaculum_thermopropionicum_SI",
        domain: "B",
        gc_percent: 52.96,
    },
    MetagenomicBin {
        id: 34,
        name: "Prochlorococcus_marinus_MIT_9313",
        domain: "B",
        gc_percent: 50.74,
    },
    MetagenomicBin {
        id: 35,
        name: "Pyrobaculum_aerophilum_IM2",
        domain: "A",
        gc_percent: 51.36,
    },
    MetagenomicBin {
        id: 36,
        name: "Ralstonia_solanacearum_PSI07",
        domain: "B",
        gc_percent: 66.13,
    },
    MetagenomicBin {
        id: 37,
        name: "Rhizobium_NGR234",
        domain: "B",
        gc_percent: 58.49,
    },
    MetagenomicBin {
        id: 38,
        name: "Rhodococcus_jostii_RHA1",
        domain: "B",
        gc_percent: 65.05,
    },
    MetagenomicBin {
        id: 39,
        name: "Rickettsia_conorii_Malish_7",
        domain: "B",
        gc_percent: 32.44,
    },
    MetagenomicBin {
        id: 40,
        name: "Rothia_dentocariosa_ATCC_17931",
        domain: "B",
        gc_percent: 53.69,
    },
    MetagenomicBin {
        id: 41,
        name: "Shigella_dysenteriae_Sd197",
        domain: "B",
        gc_percent: 51.25,
    },
    MetagenomicBin {
        id: 42,
        name: "Synechococcus_CC9605",
        domain: "B",
        gc_percent: 59.22,
    },
    MetagenomicBin {
        id: 43,
        name: "Synechococcus_JA_2_3B_a_2_13_",
        domain: "B",
        gc_percent: 58.45,
    },
    MetagenomicBin {
        id: 44,
        name: "Thermoplasma_volcanium_GSS1",
        domain: "A",
        gc_percent: 39.92,
    },
    MetagenomicBin {
        id: 45,
        name: "Treponema_pallidum_Nichols",
        domain: "B",
        gc_percent: 52.77,
    },
    MetagenomicBin {
        id: 46,
        name: "Tropheryma_whipplei_TW08_27",
        domain: "B",
        gc_percent: 46.31,
    },
    MetagenomicBin {
        id: 47,
        name: "Xenorhabdus_nematophila_ATCC_19061",
        domain: "B",
        gc_percent: 44.15,
    },
    MetagenomicBin {
        id: 48,
        name: "Xylella_fastidiosa_Temecula1",
        domain: "B",
        gc_percent: 51.78,
    },
    MetagenomicBin {
        id: 49,
        name: "_Nostoc_azollae__0708",
        domain: "B",
        gc_percent: 38.45,
    },
];

/// Preset training payloads for the metagenomic bins.
///
/// This submodule exposes get_preset_training(id) -> Option<&'static Training>,
/// where Training matches the C initialize_metagenome_X structs. We'll fill
/// them incrementally, starting with a few common bins; the rest return None
/// until ported.
pub mod presets;

pub use presets::get_preset_training_ref;
