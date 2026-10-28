fig3_fill_sensitivity.csv: calls of each Figure 3 set that are geNomad viral with a Metabolism KO, the dark-green fill, under the primary rule and five sensitivity variants.

Sets: the nine tiers of panel A after the primary-tier filter (panel 'A, primary-tier filter'), the three focal tiers of panel B in full (panel 'B focal, full tier') and the nine matched levels in full (panel 'matched level, full tier'), 21 sets and 126 rows. Tiers are cumulative, so each tier includes the calls of the stricter tiers of the same tool, and DRAM-V (aux 1), (aux 2) and (aux 3) hold the calls of the default DRAM-V configuration at auxiliary score 1, at most 2 and at most 3, so DRAM-V (aux 3) is DRAM-V (default).

primary: the fill drawn in the figure. A geNomad-viral call whose shared KofamScan KO, the best hit at e-value <= 1e-5, sits anywhere under KEGG Metabolism.

nucleotide excluded, strict: as primary, but a KO that KEGG places under Nucleotide metabolism is removed even when it is also placed under another Metabolism category.

nucleotide excluded, lenient: as primary, but a KO is removed only when Nucleotide metabolism is its only Metabolism category.

cautious excluded, CheckAMG list on shared KO: as primary, but a KO that the CheckAMG cautious-function list (AMG_filters.tsv) flags as essential, glucan, lipid, methyl or nucleotide is removed.

cautious excluded, calling tool's own flag: as primary, but a call is removed when the annotation the calling tool itself used is flagged cautious.

KofamScan adaptive threshold: a geNomad-viral call with any KofamScan hit above the adaptive threshold whose KO sits under KEGG Metabolism.

Shares are of all calls in the set. The variants were defined and compared after the primary counts had been seen, so they show how much the counts depend on the class definition and do not replace the primary rule.
