# Topics-Based Output Implementation

This directory contains documentation for the Nextflow topics-based output implementation in nf-core/rnaseq.

## Quick Start

### What are Nextflow Topics?

Nextflow 25.10+ introduces **topics** - a pub/sub mechanism that allows processes to publish outputs to named channels that can be subscribed to anywhere in the workflow, bypassing traditional channel wiring through subworkflow layers.

### Example

**Before (Traditional):**
```groovy
// In process
output:
path "*.bam", emit: bam

// In subworkflow A
emit:
bam = PROCESS.out.bam

// In subworkflow B
emit:
bam = SUBWORKFLOW_A.out.bam

// In main workflow
emit:
bam = SUBWORKFLOW_B.out.bam

// In main.nf
output:
'results' {
    WORKFLOW.out.bam >> { ... }
}
```

**After (With Topics):**
```groovy
// In process
output:
path "*.bam", emit: bam, topic: 'alignment/bam'

// In main.nf output block (directly!)
output:
'results' {
    Channel.topic('alignment/bam') >> { ... }
}
```

## Documentation Files

- **[topics_schema.md](topics_schema.md)** - Complete topic taxonomy and naming conventions
- **[topics_implementation_summary.md](topics_implementation_summary.md)** - Detailed implementation guide
- **This file** - Quick reference

## What Was Implemented

### 1. Unified Versions Collection

All 58 modules publish versions to a single `versions` topic, collected once at the end.

### 2. Direct Output Publishing

Processes publish directly to topics subscribed to in the main.nf output block.

### 3. Topic Taxonomy

Hierarchical naming: `{tool}/{datatype}/{subtype}`

Examples: `star/bam/genome`, `hisat2/logs/summary`, `umi/dedup/bam`

## Benefits

- ✅ Simplified code - no version channel mixing
- ✅ Self-documenting topic names
- ✅ Easier maintenance
- ✅ Backwards compatible

## Files Modified

- **Documentation**: 3 new files in docs/
- **Workflow**: main.nf (added output block)
- **Modules**: 58 modules updated

## References

- [Nextflow Topics Documentation](https://www.nextflow.io/docs/latest/reference/channel.html#topic)
- [Workflow Output Syntax](https://www.nextflow.io/docs/latest/workflow.html#workflow-output)
- [nf-core/rnaseq PR #1679](https://github.com/nf-core/rnaseq/pull/1679)
