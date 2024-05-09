# pollux
Optimalization of pipeline for variant calling fo respiratory samples


## scripts
`merge_samples_and_results.ipynb` - Zlepenie vsetkych google sheetov 'Respiracne vysledky' a vytvorenie z nich `respiracne_vysledky_all.tsv`. Dalej pre kazdy virus filtrovanie podla coverage, ktory vyuzije hna analyzy, ale nakoniec vyuzivame len vzorky z ktorych boli vytvarane MBS consensus sekvencie s ktorymi to budeme porovnavat.
`aligned_consensus_counts.py` - zarovanie sekvencii podla virusu, vytvorenie tabulky pre danu vzorku zo skore zarovnania (#match, #mismatch, #gap, #N, skore zarovania, skore aggregovane = skore - #gap - #mismatch pri vyuziti pariwwise2 aligneru)
`align_consensus_score_metrics.py` - ziskanie metrik zo zarovani sekvencii (m_baseline consensus vs sample consensus)

`align_cosensus_counts.py` - vytvorenie metrik zo zarovani (pairwise2.align.globalxx) -->  ._score_metrics.json
-- pre h3n2 a h1n1, lebo obsahiujeu kratke fragmety a vieme to teda 'postupne (po fragmentoch)' zarovnavat, na zarovnanie vyuzity Bio.pairwise2.align
-- nemozno pouzit pre rsv_A a rsv_B (lebo 13000nt dlhy jeden konsenzuz, tak to padalo), preto vytvoreny `aligner_score_comparison.py`, v ktorom sa vyuziva Bio.Align.PairwiseAligner, ktory zvlada aj dlhe konsenzs sekvencie rsv_a rsv_b

`aligned_PairwiseAligner_consensus_counts.ipynb` - vytvorenie metrik zo zarovani (pairwiseAligner) 

`create_configs.py` - generovanie yaml configov podla zmeny parametrov

-- spustanie congigov v screee: kazde zvlast v screene, alebo postupne v screene (prip postupne)

zostavenie tabulky: pridat zakazdym vysledky metrik (z `align_consensus_score_metrics.py`) na koniec tabulky

`create_link_to_samples.py` - create symbolic links to samples. Just sample names are changed to be later used as for configs changes (parameter in config changed) 

`n_count_in_consensus_baseline.py` - count N 'bases' in MBS consensus sequences

`aligner_score_comparison.py` - pre rozne skore zarovania vytvor pre vsetky vzniknute consenzus sekvencie (podla nastaveni parametrov) zarvonania a z nich tabulku so skore zarovnanim (skore pre match, mismatch, gap, N_count) 
-- zarovania pomocou Bio.Align.PairwiseAligner

`results_params_concat_viruses.ipynb` - zlepenie tabuliek s vysledami (vystup z aligner_score_comparison.py)

`alignment_score_plots.ipynb` - generovanie grafov vizualizujucich skore po zarovnani konsezus sekvencii s MBS konsezus sekvenciami
    - violin ploty pre virusy zvlast, violin pre virusy zvlast ale v vsetky v jednom pdf
    - swarm plot pre vsetky virusy, aby sme ukazali, ze toto nastavenie je dobre u vsetkych virusoch a ako sa meni skore ak pouzijeme iny parameter

`alignment_score_plots_all_viruses.ipynb` 
- generovanie grafu pre vsetky virusy (skore zarovnania a jednotlive skumane parametre)  -tieto zlepit potom do jedneho grafu
- dalsia vizualizacia porovnania skore matica(virus, paramtere): med(v1, q=30)/med(v1, q=MBS); tiez potom zlepit tieto grafy

`merge_aligned_scores_all_viruses.ipynb` - zlepenie vystupov generovanych skriptom `aligner_score_comparison.py` : to je pre kazdy virus zvlast subor, preto tu je pre danu skorovaciu schemu zlepene skore pre vseky virusy (do jedneho suboru)
- pridanie stlpca oznacujuce ten riadok, ktore predstavuje to vyhodnotenie, ktore predstavuje MBS (porovnavane samo so sebou)

`unique_variant_sars-cov2_nexstrainclades.py` - get unique variants according to nexstrain clades classification system from LAPIS

`unique_variant_sars-cov2_pangolineage.py` - get unique variants (lineage) according to pango  classification system from LAPIS - run it on server (gen-manager)

`unique_variant_sars-cov2_pangolineage_old.py` - get unique variants (lineage) according to pango  classification system from LAPIS

`aligner_score_comparison_cutadapt.py` - to iste ako aligner_score_comparison.py ale pre trimmer nastroj cutadapt



`create_configs_all_samples_cutadapt_mbs.py` - nastaveie mbs pre cutadapt, aj mbs nastavenia pre variant callery: generuje konfiguracne skripty v cutadapt_mbs_vc_mbs_configs/
-- vysledne zarovannia v subore: `results_all_viruses_match1_mis-1_o-2_e-2_2024_4_15_cutadapt_mbs_vc_mbs.tsv`

`create_list_parameter_set.ipynb` - vytvorenie df pre zoznam testovanych parametrov (takym nazvom su oznacene vysledne vzorky), ktore su dalej pouzivane v  `aligner_score_comparison....py`

Generovanie mbs nastaveni cutadaptu a mbs nastaveni variant callerov:
`create_configs_all_samples_cutadapt_mbs.py` generoval pomocou konfigu `[nazov variant calleru]_config_cutadapt_mbs.yaml` konfigy ulozene v : `scripts/variant_call/cutadapt_mbs_vc_mbs_configs/...yaml`
`scripts/create_links/create_link_to_all_samples_cutadapt_mbs_vc_mbs.py`
`aligner_score_comparison_cutadapt_mbs_vc_mbs.py` - tieto vysledky zmergovane v `merge_aligned_scores_all_virus.ipynb` a vytvoreny subor `./data/results_metrics/results_all_viruses/results_all_viruses_match1_mis-1_o-2_e-2_2024_4_15_cutadapt_mbs_vc_mbs.tsv`

`aligner_score_comparison_variantcall.py`(ako aj `.ipynb`) - zarovanie s baseline consenzus sekvenciou, a variant callermi, pri ktorych sa menili testovane paramtre
(configy pre beh analyz vytvorene v `create_configs_[nazov variant callera]_variantc_all_samples.py`)

`aligner_score_comparison_trimm_vc.py` (ako aj `.ipynb`) - zarovanie s baseline consenzus sekvenciou, trimm = trimmers(trimmomatic a cutadapt), VC = variant_callers, pre obe mnoziny nastrojov vybrane naj nastavenie parametrov a ich kombinacia
(configy pre beh analyz vytvorene v `create_configs_all_samples_final_comparison.py`)

`merge_aligned_scores_all_virus.ipynb` 
-- tiez obsahuje zlepenie: finalne vythodnotenia (trimmer + variant callers oba s naj nastavenim) so skore s mbs nastavenim pre dany variant caller a mbs nastavebim trimmeru ==> vystup: `results_all_viruses_match1_mis-1_o-2_e-2_2024_4_15_cutadapt_vc_mbs.tsv` a `results_all_viruses_match1_mis-1_o-2_e-2_2024_4_14_trimmomatic_vc_mbs.tsv`, z ktorych su generovane heatmapy v ktorych je porovanie ako sa zmenilo skore zarovanania z povodneho mbs nastavenia trimmeru a mbs nastavenia VC - zmena skore pri pouziti vybranych naj nastaveni trimmeru a VC 

Pre spracovanie a vygenerovanie grafov spusti:
1. `python aligner_score_comparison_cutadapt.py` alebo `aligner_score_comparison.py`
2. `merge_aligned_scores_all_viruses.ipynb` - zlep vystupy
3. `alignment_score_plots_all_viruses.ipynb` - vygeneruj grafy

## directory structructure

`scripts_main.zip` -- generovanie skriptov a vysledkov

`scripts.zip` -- vygenerovane skripty (konfiguracne subory, skripty pre spustanie konfiguracnych suborov)

`DP-respiracne.zip` -- baseline konsezus sekvencie, a data (vysledne tabulky, vyhodnotenia zarovani, z vysledkov vygenerovane grafy) 


`consensus_baseline` - konsenzus sekvencie (skontrolovane Mirom)
`data`
    - `vysledky` - tsv z google sheetov (Respiracne_vysledky)
    - `report_pollux` - namountovany priecinok /data/project/pollux/report/public/ z gen-manageru 
    - `result_metrics` - vyhodnotenia zarovnani (vegenerovane `aligned_consensus_counts.py`)

`scripts/variant_call` - configs for qsnake (ivar.yaml is baseline config). Those directories have to be created before running `create_configs.py` script.
    - `rsvA_configs` - configs with changed parameters. Generated from `create_configs.py`
    - `h1n1_configs` 


