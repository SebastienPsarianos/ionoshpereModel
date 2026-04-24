# Thesis Revision TODO List

## Chapter 3: Numerical Implementation

### §3.1 Discretization


- [ ] "`% TODO: Determine the normal vector difference and include it in this equation`" comment after Eq. `eq:pole_integration_exp` — resolve or delete

### §3.2.1 Discretization of Polar BC

- [ ] **Periodicity in $\phi$ not mentioned anywhere** — add one sentence stating that indices $(a, i \pm 1)$ wrap modulo $n_\phi$ at $i = 0$ and $i = n_\phi - 1$
- [ ] "first-derivative central-difference approximation" — confirm hyphenation consistency throughout section (some instances of "first derivative" without hyphen still appear)
- [ ] The `\lim_{\theta\to 0}` argument is technically showing that $\partial_\phi\psi = O(\theta)$; currently framed as "in the limit" which is fine, but worth clarifying that $\theta_0 = \Delta\theta/2$ is small but non-zero in the actual discretization
- [ ] Fig. `fig:boundary_condition_contour` caption: "southward interplanetary magnetic field" is good; confirm IMF is spelled out at first use and abbreviated thereafter if it appears elsewhere

### §3.3 Solver

#### §3.3.1 Tpetra

- [ ] "trilinos" (lowercase) in "All modern trilinos packages" — capitalize as "Trilinos"
- [ ] **Overclaim:** "All modern trilinos packages rely on Tpetra as their common interface layer and therefore it is a requirement for any Trilinos solve." — soften to something like "Most modern Trilinos solver and preconditioner packages, including Belos and MueLu used in this work, operate on Tpetra objects."
- [ ] "This allows the ability to tweak the distribution on a global scale" — "allows the ability to" is redundant; simplify to "This allows the distribution to be modified globally by changing a single `Tpetra::Map` instance."
- [ ] "a sequential rows $\phi$-major map" — ambiguous phrasing; clarify as "a $\phi$-major ordering with contiguous blocks of rows assigned to each MPI rank"
- [ ] `Zoltan~(\cite{ZoltanUsersGuide})` — extra parentheses around `\cite`; should be `Zoltan~\cite{ZoltanPaper}` (and consider switching to the Boman et al. 2012 *Scientific Programming* paper rather than the user's guide)
- [ ] "interleaved storage between ranks" — "interleaved distributions" reads more naturally than "interleaved storage"
- [ ] Paragraph break issue: `\\` at end of first paragraph of §3.3.1 is non-standard for paragraph separation; use a blank line instead

#### §3.3.2 Belos

- [ ] "belos solver" (lowercase) in "The belos solver and the associated solver parameters" — capitalize as "Belos"
- [ ] "`\texttt{GMRES}` solver" — GMRES is an algorithm name, not a code identifier; use plain text "GMRES solver" unless referring specifically to the SolverFactory string
- [ ] **Strengthen the GMRES justification:** "this showed the highest performance on the nonsymmetric matrix" — the Belos paper frames GMRES as favored for *robustness* over short-recurrence methods on nonsymmetric systems; consider: "GMRES with preconditioning was selected as it showed the best convergence behaviour for the nonsymmetric linear system; Belos favours GMRES variants for nonsymmetric systems on robustness grounds over short-recurrence methods such as TFQMR~\cite{BelosPaper}."
- [ ] **Add explanation of WHY the system is nonsymmetric** — one sentence attributing it to the first-derivative terms and the Hall conductance cross-coupling in Eq. `eq:expanded_continuity_equation`
- [ ] Consider mentioning nonsymmetry earlier (end of §3.1 or §3.2) with a forward-reference to §3.3.2, so the solver choice is motivated when the reader first encounters it

#### §3.3.3 Preconditioning

- [ ] **"IfPack" — likely should be "Ifpack2"** (the modern Tpetra-compatible version). "Ifpack" (without the 2) refers to the older Epetra version. Verify against your actual build (`spack.yaml` or source code grep).
- [ ] "incomplete lower upper triangularization" → "incomplete LU (ILU) factorization" (standard terminology)
- [ ] "were found to be optimal through empirical testing" → "were found to be most effective" or "best-performing" (can't prove optimality from finite empirical testing)
- [ ] Consider adding an Ifpack2 citation if you confirm you're using Ifpack2

#### §3.3.4 Application Structure

- [ ] "discuessed" → "discussed" (typo)
- [ ] "recieving" → "receiving" (typo: i before e except after c)
- [ ] "feasability" → "feasibility" (typo)
- [ ] "the figures in this paper" → "the figures in this thesis" (or "chapter")
- [ ] "ionosphere model used in~\cite{global_groth_2000}" → "ionosphere model described in~\cite{global_groth_2000}"
- [ ] The ghost-ring mention returns here as "in the process of implementation" — confirm this wording doesn't contradict the Polar BC section framing
- [ ] Fig. `diag:solver_flow` caption: "configurable YAML parameters" — drop "YAML" since it's already clear from context and the diagram itself
- [ ] `\\` at paragraph end is non-standard for paragraph separation; use blank lines
- [ ] Paragraph starting "The MHD interpolation is currently in progress" — "circular loop structure" is tautological; "loop structure" alone suffices
- [ ] Filename `ApplicationFLow` (capital L) — likely typo for `ApplicationFlow`; confirm the file in `diagrams/` matches

### §3.4 Convergence Study and Solver Optimization

- [ ] **All figure captions are placeholders** ("TEST1", "TEST", "TEST2", empty captions) — populate these
- [ ] **Section has no introductory paragraph** — add a lead-in explaining what's being studied, across what parameters, using what problem size / configuration
- [ ] Fig. `fig:overhead` caption ends mid-sentence: "Mean radial current processing and overhead for the MueLu" — finish it
- [ ] Fig. `diag:full_flow` is currently at the end of §3.4, structurally misplaced — belongs near §3.3.4 or moved to appendix (see below)
- [ ] Fig. `diag:full_flow` caption ends mid-sentence: "configurable parameters with" — finish it
- [ ] Fig. `diag:full_flow` caption starts identically to Fig. `diag:solver_flow` ("Application structure diagram displaying the data flow...") — differentiate them
- [ ] Filename `ApplicationFLow` typo (as above)
- [ ] Fig. `fig:time_breakdown` and Fig. `fig:scaling_wallclock` have empty captions — fill in
- [ ] Consider moving `diag:full_flow` to an appendix with forward-references from the intro and Chapter 3 opening

## Global / Cross-cutting issues

### Terminology consistency

- [ ] "Pederson" vs "Pedersen" — audit all instances, use "Pedersen"
- [ ] "first derivative" vs "first-derivative" (as compound modifier) — audit hyphenation
- [ ] "five point stencil" vs "five-point stencil" — hyphenate when used as a compound modifier
- [ ] "finite-volume style" vs "finite-volume-style" — hyphenate fully when compound modifier
- [ ] "Trilinos" capitalization — proper noun, always capitalized
- [ ] "Belos", "Tpetra", "MueLu", "Ifpack2", "Kokkos" — all proper package names, always capitalized as shown
- [ ] "magnetohydrodynamics" / "MHD" — spell out on first use in each chapter if not done earlier
- [ ] "IMF" — spell out on first use ("interplanetary magnetic field")
- [ ] "FAC" — spell out on first use ("field-aligned current")
- [ ] "SZA" and "MLT" — spell out on first use if they appear in prose (not just diagrams)

### Notation consistency

- [ ] $\psi$ (potential) — confirm used consistently throughout
- [ ] $\Sigma$ vs $\sigma$ for conductances — $\Sigma$ appears to be the convention; audit for stray $\sigma$
- [ ] $R_E$ — confirm consistent subscript (not $R_e$); Eq. `eq:pole_integration_line2` uses $R_e$ (lowercase e), should be $R_E$
- [ ] $\theta_0$ definition — used for both north ($\Delta\theta/2$) and south ($\pi - \Delta\theta/2$) pole; consider adding a subscript or note for clarity
- [ ] Superscript conventions $(i,j)$, $(a,i)$, $(\parC, i)$, $(\rm pole)$ — all introduced in §3.2.1 but worth confirming they're used consistently throughout

### Citation audit

- [ ] Verify every `\cite` resolves to an entry in the `.bib` file (common for citations to break during revision)
- [ ] `ZoltanUsersGuide` — consider replacing with `ZoltanPaper` (Boman et al. 2012)
- [ ] Confirm Ifpack2 citation added if using Ifpack2
- [ ] `TrilinosPaper`, `TpetraPaper`, `BelosPaper`, `MueLuUsersGuide` — all present and resolving
- [ ] `global_groth_2000` — confirm exists
- [ ] `goodman1995`, `amm1996`, `William_H_Press1992-ps` — confirm all present

### Figure / label issues

- [ ] `diag:full_flow` caption incomplete (sentence ends at "with")
- [ ] `fig:overhead` caption incomplete (sentence ends at "MueLu")
- [ ] `fig:scaling_time` caption is "TEST1" (placeholder)
- [ ] `fig:scaling_iters` caption is "TEST" (placeholder)
- [ ] `fig:scaling_time_iters` caption is "TEST2" (placeholder)
- [ ] `fig:scaling_wallclock` caption empty
- [ ] `fig:time_breakdown` caption empty
- [ ] All `../convergenceAnalysis/plotsAndData/` image paths — confirm these resolve correctly from the thesis build directory
- [ ] Figure placement: `diag:full_flow` misplaced in §3.4
- [ ] Consider if Fig. `fig:boundary_condition_contour` path `figures/boundary condition zoomed.png` should have underscores instead of spaces (spaces in paths can cause LaTeX headaches)

### Structural

- [ ] Decide final location for `diag:full_flow`: appendix (recommended) vs Chapter 3 intro vs keep in place
- [ ] Add forward-reference from intro to appendix if moving diagram there
- [ ] Consider whether §3.3 needs an opening paragraph introducing its three subsubsections before diving into Tpetra
- [ ] Confirm §3.3 subsubsection nesting (`\subsubsection{Tpetra}` etc.) matches your chapter's sectioning conventions

### Writing style (revision-phase sweep)

- [ ] "allows the ability to" pattern — redundant, simplify throughout
- [ ] "through defining", "through the use of" — often simplifiable to "by defining", "by using"
- [ ] Passive voice audit — not required to eliminate, but some sentences gain clarity from active voice
- [ ] Paragraph-separator `\\` usage — should generally be blank lines for paragraph breaks; `\\` is for line breaks within a paragraph
- [ ] Read aloud (or use TTS) during final revision — catches typos and awkward phrasings the eye skips

### Pre-submission checks

- [ ] Resolve all `% TODO:` comments
- [ ] Full citation audit (every `\cite` and every `\ref` resolves)
- [ ] Compile check with no errors or warnings
- [ ] Page-break audit (no orphan headers, no awkward figure placements)
- [ ] Final read-through as a reader, not a writer
