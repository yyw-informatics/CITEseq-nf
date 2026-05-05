/*
 * Branch Selection: Promote decontamination winner to canonical paths
 *
 * Contract:
 *   Input:  Branch name (e.g. "decontx") + completion signal
 *   Action: Copies adata_qc_{branch}.h5ad -> adata_qc.h5ad (and adt_raw)
 *   Output: Canonical adata_qc.h5ad and adata_adt_raw.h5ad paths
 *
 * Runs after both decontamination branches complete step01b (or
 * after a single branch if only one is used).
 *
 * Channel-pure: no publishDir. Publishing is done at workflow level.
 */

process select_branch {
    label 'cpu_small'
    label 'env_cite_seq'

    input:
    val(branch)
    val(ready)

    output:
    val(true), emit: done

    script:
    """
    ${params.python} ${workflow.projectDir}/scripts/admin/select_branch.py ${branch}
    """
}
