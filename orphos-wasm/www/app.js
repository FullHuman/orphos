import init, { init_panic_hook, analyze_sequence } from './pkg/prodigal_wasm.js';

// Initialize WASM module
let wasmInitialized = false;

async function initWasm() {
    if (!wasmInitialized) {
        await init();
        init_panic_hook();
        wasmInitialized = true;
        console.log('WASM module initialized');
    }
}

// DOM Elements
const fastaFileInput = document.getElementById('fasta-file');
const fileNameDisplay = document.getElementById('file-name');
const fastaTextarea = document.getElementById('fasta-input');
const analyzeBtn = document.getElementById('analyze-btn');
const loadingDiv = document.getElementById('loading');
const resultsSection = document.getElementById('results-section');
const errorSection = document.getElementById('error-section');
const errorMessage = document.getElementById('error-message');
const outputPre = document.getElementById('output');
const geneCountSpan = document.getElementById('gene-count');
const sequenceCountSpan = document.getElementById('sequence-count');
const downloadBtn = document.getElementById('download-btn');
const copyBtn = document.getElementById('copy-btn');

// Options
const modeSelect = document.getElementById('mode');
const formatSelect = document.getElementById('format');
const translationTableSelect = document.getElementById('translation-table');
const closedEndsCheckbox = document.getElementById('closed-ends');
const maskNRunsCheckbox = document.getElementById('mask-n-runs');
const forceNonSdCheckbox = document.getElementById('force-non-sd');

// File input handler
fastaFileInput.addEventListener('change', async (e) => {
    const file = e.target.files[0];
    if (file) {
        fileNameDisplay.textContent = file.name;
        try {
            const text = await file.text();
            fastaTextarea.value = text;
        } catch (error) {
            showError(`Failed to read file: ${error.message}`);
        }
    } else {
        fileNameDisplay.textContent = 'No file chosen';
    }
});

// Analyze button handler
analyzeBtn.addEventListener('click', async () => {
    const fastaContent = fastaTextarea.value.trim();
    
    if (!fastaContent) {
        showError('Please provide a FASTA sequence');
        return;
    }

    if (!fastaContent.startsWith('>')) {
        showError('Invalid FASTA format. Sequence must start with ">"');
        return;
    }

    try {
        // Initialize WASM if not already done
        await initWasm();

        // Show loading state
        showLoading(true);
        hideError();
        hideResults();

        // Get options
        const options = {
            mode: modeSelect.value,
            format: formatSelect.value,
            closed_ends: closedEndsCheckbox.checked,
            mask_n_runs: maskNRunsCheckbox.checked,
            force_non_sd: forceNonSdCheckbox.checked,
            translation_table: translationTableSelect.value ? parseInt(translationTableSelect.value) : null
        };

        console.log('Analysis options:', options);

        // Run analysis
        const startTime = performance.now();
        const result = analyze_sequence(fastaContent, options);
        const endTime = performance.now();

        console.log(`Analysis completed in ${(endTime - startTime).toFixed(2)}ms`);

        // Display results
        displayResults(result);

    } catch (error) {
        console.error('Analysis error:', error);
        showError(error.message || error.toString());
    } finally {
        showLoading(false);
    }
});

// Download button handler
downloadBtn.addEventListener('click', () => {
    const content = outputPre.textContent;
    const format = formatSelect.value;
    const extension = format === 'gbk' ? 'gbk' : 
                     format === 'gff' ? 'gff' : 
                     format === 'sco' ? 'sco' : 'gca';
    
    const blob = new Blob([content], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `prodigal_results.${extension}`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
});

// Copy button handler
copyBtn.addEventListener('click', async () => {
    const content = outputPre.textContent;
    try {
        await navigator.clipboard.writeText(content);
        const originalText = copyBtn.textContent;
        copyBtn.textContent = '✓ Copied!';
        setTimeout(() => {
            copyBtn.textContent = originalText;
        }, 2000);
    } catch (error) {
        showError('Failed to copy to clipboard');
    }
});

// Helper functions
function showLoading(show) {
    loadingDiv.style.display = show ? 'flex' : 'none';
    analyzeBtn.disabled = show;
}

function showError(message) {
    errorMessage.textContent = message;
    errorSection.style.display = 'block';
}

function hideError() {
    errorSection.style.display = 'none';
}

function hideResults() {
    resultsSection.style.display = 'none';
}

function displayResults(result) {
    outputPre.querySelector('code').textContent = result.output;
    geneCountSpan.textContent = result.gene_count;
    sequenceCountSpan.textContent = result.sequence_count;
    resultsSection.style.display = 'block';
    
    // Scroll to results
    resultsSection.scrollIntoView({ behavior: 'smooth' });
}

// Load example sequence on page load (optional)
window.addEventListener('DOMContentLoaded', () => {
    const exampleSequence = `>Example_sequence
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCG
GGCTGAGATCGCAATGAAGAAGGTGGCTTTTGCTACTGACGGTGATTACCCGAGTTATCA
CCAATATCGTGTACAACCGTTGAATTATGCGATCATCGCCGCCGATGCCACCATGGCCGC
TGAAGCGGCTGGCGCCGACGCTGTCCGTCGCCTGATCCAGGCCTTTATCGCCGGTATTTC
GCACGAGCTTCCGCATGCCGAAATCGCCCAGCGCGGTGAGTTACAGGCCTATCGTGGTCA
TCCGAAAACTTGA`;
    
    // Uncomment to load example by default
    // fastaTextarea.value = exampleSequence;
});

// Initialize WASM on page load
initWasm().catch(err => {
    console.error('Failed to initialize WASM:', err);
    showError('Failed to initialize WebAssembly module. Please refresh the page.');
});
