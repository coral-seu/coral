# coral

## User‘s Guide

This script is designed for clustering coral single-cell transcriptome data. After researchers obtain the coral single-cell transcriptome results and the expression matrix information, they will first perform initial processing using the Seurat package:

1. Initially, CCA (Canonical Correlation Analysis) is used for dimensionality reduction, followed by L2 normalization. The Most Neighbors (MNN) method is then employed to find anchors (similar cell pairs), anchor scores (SNN) are calculated, and low-quality anchors are filtered out to integrate the data.
2. The FindVariableGenes function is used to identify variable gene information and find the genes with the highest scores.
3. The ScaleData function is applied to scale the data, normalizing it for further analysis.
4. The RunPCA function is utilized to perform principal component analysis (PCA) on the processed data, mapping high-dimensional cell representation information into a lower-dimensional space.
5. The JackStraw function is employed to determine the statistical significance of the PCA scores.
6. The K-Nearest Neighbor (KNN) algorithm is used to calculate a network graph with cells as nodes. Then, based on the Louvain community detection algorithm, the cells are divided into corresponding clusters.
7. The results are visualized using t-SNE or UMAP methods, preserving the global structure of the data.

Potential marker gene information is obtained based on the highly expressed genes.

Subsequently, using the marker gene information provided by this script, researchers can manually match and identify the cell types corresponding to the clustered clusters.

## Manual annotation of cell type

By combining the data obtained from the experimental procedures and referencing relevant literature sources, this script offers a series of marker gene information, which typically shows high expression in specific cell types of coral or other related species. At the end of each gene, a notation is added to indicate the source of the marker gene information, whether it is derived from direct data from this experiment, published literature, or other reliable database resources. 

| Cell  Type                 | Cluster number   | Symbol           | Gene id                | Marker gene                                                  | Orthologous gene                                             | Reference                                                    |
| -------------------------- | ---------------- | ---------------- | ---------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 1. Cnidocytes              | C5, C9, C12, C23 | LRX3             | evm.model.Chr12.438    | PREDICTED: leucine-rich repeat extensin-like protein 3 [Acropora digitifera] | XP_015764414.1                                               | This project                                                 |
|                            |                  | LRX6             | evm.model.Chr12.437    | leucine-rich repeat extensin-like protein 6 [Acropora millepora] | XP_029208639.1                                               | This project                                                 |
|                            |                  | COL-2            | evm.model.Chr14.183    | cuticle collagen 2C-like   [Acropora millepora]              | XP_044169233.1                                               | This project                                                 |
|                            |                  | CPD              | evm.model.Chr5.680     | carboxypeptidase D-like   [Acropora millepora]               | XP_029191340.2                                               | This project                                                 |
|                            |                  | Minicollagen_1   | evm.model.Chr12.542    | minicollagen                                                 | Nv-Ncol-5, Ms-Ncol-4                                         | https://doi.org/10.1016/j.tig.2008.07.001                    |
|                            |                  | Minicollagen_2   | evm.model.Chr12.543    | minicollagen                                                 | Nv-Ncol-5, Ms-Ncol-4                                         | https://doi.org/10.1016/j.tig.2008.07.001                    |
|                            |                  | Minicollagen_3   | evm.model.Chr12.357    | minicollagen                                                 | Ad-Ncol-1                                                    | https://doi.org/10.1016/j.tig.2008.07.001                    |
|                            |                  | Nematogalectin_1 | evm.model.Chr9.986     | nematogalectin                                               | Hm-NemgalA, Ha-NemgalA, Ho-NemgalA, Hv-NemgalA, Hm-NemgalB,  Ha-NemgalB, Ho-NemgalB, Ch-NemgalA, Ch-NemgalB, Au-Nemgal, Ms-Nemgal,  Ap-Nemgal, Av-Nemgal, Nv-Nemgal, Ad-Nemgal, Xe_006067 | https://doi.org/10.1073/pnas.1003256107     https://doi.org/10.1038/s41586-020-2385-7 |
|                            |                  | Nematogalectin_2 | evm.model.Chr9.985     | nematogalectin                                               | Xe_006066, v1g232014, v1g232015                              | https://doi.org/10.1038/s41586-020-2385-7     https://doi.org/10.1016/j.cell.2018.05.019 |
|                            |                  | Nematogalectin_3 | evm.model.Chr9.982     | nematogalectin                                               | Hm-NemgalA, Ha-NemgalA, Ho-NemgalA, Hv-NemgalA, Hm-NemgalB,  Ha-NemgalB, Ho-NemgalB, Ch-NemgalA, Ch-NemgalB, Au-Nemgal, Ms-Nemgal,  Ap-Nemgal, Av-Nemgal, Nv-Nemgal, Ad-Nemgal, Xe_006067 | https://doi.org/10.1073/pnas.1003256107     https://doi.org/10.1038/s41586-020-2385-7 |
|                            |                  | Nematogalectin_4 | evm.model.Chr9.1084    | nematogalectin                                               | Hm-NemgalA, Ha-NemgalA, Ho-NemgalA, Hv-NemgalA, Hm-NemgalB,  Ha-NemgalB, Ho-NemgalB, Ch-NemgalA, Ch-NemgalB, Au-Nemgal, Ms-Nemgal,  Ap-Nemgal, Av-Nemgal, Nv-Nemgal, Ad-Nemgal, Xe_006067 | https://doi.org/10.1073/pnas.1003256107     https://doi.org/10.1038/s41586-020-2385-7 |
|                            |                  | POU4F3           | evm.model.Chr5.322     | POU domain, class 4, transcription factor 3                  | v1g160868                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
| 2.  Neurons 1              | C0, C10, C11     | AKAP9            | evm.model.Chr10.85     | A-kinase anchor protein 9-like isoform X1 [Acropora millepora] | XP_029186357.2                                               | This project                                                 |
|                            |                  | GCLM_1           | evm.model.Chr10.1143.1 | Glutamate--cysteine ligase regulatory subunit [Acropora  cervicornis] | KAK2555327.1                                                 | This project                                                 |
|                            |                  | GCLM_2           | evm.model.Chr7.11      | glutamate--cysteine ligase catalytic subunit-like [Acropora  millepora] | XP_029207154.2                                               | This project                                                 |
|                            |                  | PKD1             | evm.model.Chr6.549     | polycystic kidney disease protein 1                          | v1g21107                                                     | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | PKD2             | evm.model.Chr5.151     | polycystic kidney disease protein 2                          | v1g160849                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | SCN1             | evm.model.Chr2.216     | Sodium channel protein 1 brain                               | NVE7349                                                      | https://doi.org/10.1016/j.cell.2018.05.019                   |
| 3. Neurons 2               | C2               | CD151            | evm.model.Chr6.1145    | CD151 antigen                                                | v1g233307                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | SYN2             | evm.model.Chr5.165     | synapsin-2                                                   | v1g110232                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | KCNC1            | evm.model.Chr6.1071    | potassium voltage-gated channel subfamily C member 1         | v1g21646                                                     | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | GATA2            | evm.model.Chr5.826     | GATA-binding factor 2                                        | XP_022781270.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | ETV1             | evm.model.Chr12.187    | ETS translocation variant 1                                  | XP_022794506.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | Rfamide          | evm.model.Chr8.494     | antho-RFamide neuropeptides                                  | v1g1374                                                      | https://doi.org/10.1016/j.cell.2018.05.019                   |
| 4. Gland cells             | C13              | ADT-1_1          | evm.model.Chr14.636    | A disintegrin and metalloproteinase with thrombospondin motifs  adt-1-like [Acropora millepora] | XP_044166601.1                                               | This project                                                 |
|                            |                  | ADT-1_2          | evm.model.Chr14.632    | A disintegrin and metalloproteinase with thrombospondin motifs  adt-1-like [Acropora millepora] | XP_044166601.1                                               | This project                                                 |
|                            |                  | MUC5AC           | evm.model.Chr11.84     | mucin-5AC-like                                               | XP_022788291.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | RFX4             | evm.model.Chr1.1331    | transcription factor RFX4-like                               | XP_022779108.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | CHST3            | evm.model.Chr14.888    | carbohydrate sulfotransferase 3-like                         | XP_022782163.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | PAPSS1           | evm.model.Chr9.540     | bifunctional 3'-phosphoadenosine 5'-phosphosulfate  synthase-like | XP_022803433.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
| 5. Epidermal cells         | C16              | MOXD1            | evm.model.Chr4.1732    | PREDICTED: DBH-like monooxygenase protein 1 homolog [Acropora digitifera] | XP_015776337.1                                               | This project                                                 |
|                            |                  | PXDN 2           | evm.model.Chr3.1094    | peroxidasin homolog                                          | XP_022798424.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | PXDN 1           | evm.model.Chr3.1095    | peroxidasin-like                                             | XP_022798424.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
| 6. Calicoblastic cells     | C4               | GOLGA4           | evm.model.Chr1.1186    | golgin subfamily A member 4-like [Acropora millepora]        | XP_015760268.1                                               | This project                                                 |
|                            |                  | MMP24            | evm.model.Chr2.747     | matrix metalloproteinase-24-like [Acropora millepora]        | XP_029201938.2                                               | This project                                                 |
|                            |                  | SLC4A8           | evm.model.Chr6.1415    | electroneutral sodium bicarbonate exchanger 1-like isoform  X1 [Acropora millepora] | XP_029209647.2                                               | This project                                                 |
|                            |                  | CADN             | evm.model.Chr13.159    | coadhesin                                                    | B3EWZ3                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst110                    |
|                            |                  | GXN2             | evm.model.Chr12.608    | galaxin-2                                                    | B8UU51                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst109                    |
|                            |                  | ASOMP            | evm.model.Chr10.1330   | acidic skeletal organic matrix protein                       | B3EWY7                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst109                    |
|                            |                  | GXN              | evm.model.Chr12.609    | galaxin                                                      | D9IQ16                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst109                    |
|                            |                  | SAAR1            | evm.model.Chr6.964     | skeletal aspartic acid-rich protein 1                        | B3EWY6                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst112                    |
|                            |                  | SAAR2            | evm.model.Chr6.888     | skeletal aspartic acid-rich protein 2                        | B3EWY8                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst111                    |
|                            |                  | TRP              | evm.model.Chr9.205     | threonine-rich protein                                       | B3EWZ7                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst112                    |
|                            |                  | Mucin            | evm.model.Chr2.1104    | Mucin-like protein                                           | B3EWY9                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst112                    |
|                            |                  | CDP              | evm.model.Chr9.207     | CUB domain-containing protein                                | B3EX01                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst112                    |
|                            |                  | Ectin            | evm.model.Chr3.309     | Ectin                                                        | B3EWZ8                                                       | https://doi.org/10.1093%2Fmolbev%2Fmst112                    |
| 7. Alga-hosting cells      | C7, C20, C22     | CA2              | evm.model.Chr3.417     | Carbonic anhydrase 2 [Acropora cervicornis]                  | KAK2574613.1                                                 | This project                                                 |
|                            |                  | SLC26A6          | evm.model.Chr12.6      | Solute carrier family 26 member 6 [Acropora cervicornis]     | KAK2573107.1                                                 | This project                                                 |
|                            |                  | DEGS1            | evm.model.Chr2.728     | Sphingolipid delta(4)-desaturase DES1 [Acropora cervicornis] | KAK2560926.1                                                 | This project                                                 |
|                            |                  | SGPL1            | evm.model.Chr9.559     | sphingosine-1-phosphate lyase 1-like [Acropora millepora]    | XP_044174562.1                                               | This project                                                 |
|                            |                  | NADP-GDH         | evm.model.Chr6.501     | NADP-specific glutamate dehydrogenase [Acropora cervicornis] | KAK2571680.1                                                 | This project                                                 |
|                            |                  | CTSL1_1          | evm.model.Chr9.1138    | cathepsin L1-like isoform X1   [Acropora millepora]          | XP_029191107.2                                               | This project                                                 |
|                            |                  | CTSB             | evm.model.Chr7.425     | cathepsin B-like   [Acropora millepora]                      | XP_022786917.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | PSAP             | evm.model.Chr4.674     | prosaposin                                                   | amur_s0134.g9.t1                                             | https://doi.org/10.1093%2Fmolbev%2Fmsaa216                   |
|                            |                  | APOD             | evm.model.Chr2.1153    | apolipoprotein D                                             | XP_022787119.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | NPC2L            | evm.model.Chr1.1138    | NPC intracellular cholesterol transporter 2-like             | /                                                            | https://doi.org/10.1111/cmi.12564                            |
|                            |                  | ALDH-NADP        | evm.model.Chr4.14      | aldehyde dehydrogenase, dimeric NADP-preferring-like         | EDO37493.1                                                   | https://doi.org/10.1111/cmi.12564                            |
| 8.  Digestive filaments    | C6               | ERG              | evm.model.Chr1.917     | Transcriptional regulator Erg [Acropora cervicornis]         | KAK2573353.1                                                 | This project                                                 |
|                            |                  | RFX6             | evm.model.Chr9.427     | DNA-binding protein RFX6-like                                | v1g85569                                                     | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | MPP7             | evm.model.Chr11.1090   | MAGUK p55 subfamily member 7                                 | v1g199299                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | Hedgehog         | evm.model.Chr12.649    | tiggy-winkle hedgehog protein                                | v1g241466                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | FoxA2A           | evm.model.Chr11.1101   | forkhead box protein A2-A-like                               | v1g165261                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
| 9. Gastrodermal cells      | C3, C8, C19      | LOC107350580     | evm.model.Chr1.873     | PREDICTED: uncharacterized protein LOC107350580 [Acropora digitifera] | XP_015772309.1                                               | This project                                                 |
|                            |                  | ELP2             | evm.model.Chr4.1609    | elongator complex protein 2-like [Acropora millepora]        | XP_044172810.1                                               | This project                                                 |
|                            |                  | SVOP             | evm.model.Chr3.463     | synaptic vesicle 2-related protein-like isoform X1 [Acropora millepora] | XP_029194443.2                                               | This project                                                 |
|                            |                  | KCNJ2            | evm.model.Chr11.853    | inward rectifier potassium channel 2-like [Acropora millepora] | XP_029204002.1                                               | This project                                                 |
|                            |                  | SAP130           | evm.model.Chr1.632     | histone deacetylase complex subunit SAP130-like isoform  X1 [Acropora millepora] | XP_029207818.2                                               | This project                                                 |
|                            |                  | COL1A2_1         | evm.model.Chr12.961    | collagen alpha-2(I) chain                                    | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A1_1         | evm.model.Chr12.964    | collagen alpha chain                                         | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A2_2         | evm.model.Chr10.1599   | collagen alpha-2(I) chain                                    | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A1_2         | evm.model.Chr10.1230   | collagen alpha-1(I) chain                                    | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A1_3         | evm.model.Chr12.947    | collagen alpha chain                                         | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A1_4         | evm.model.Chr12.736    | collagen alpha-1(I) chain                                    | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A1_5         | evm.model.Chr12.963    | collagen alpha-1(I) chain                                    | v1g80945                                                     | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A2_3         | evm.model.Chr12.957    | collagen alpha-2(I) chain                                    | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A1_6         | evm.model.Chr12.954    | collagen alpha-1(I) chain                                    | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A1_7         | evm.model.Chr12.958    | collagen alpha-1(I) chain                                    | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | COL1A1_8         | evm.model.Chr12.948    | collagen alpha chain                                         | v1g204742                                                    | https://doi.org/10.1016/j.cell.2018.05.019                   |
| 10. Retractor muscle cells | C17              | LOC114955083     | evm.model.Chr9.447     | uncharacterized protein LOC114955083 [Acropora millepora]    | XP_029187668.2                                               | This project                                                 |
|                            |                  | TIMP3            | evm.model.Chr5.183     | PREDICTED: metalloproteinase inhibitor 3-like [Acropora digitifera] | XP_015769572.1                                               | This project                                                 |
|                            |                  | MMP2             | evm.model.Chr2.746     | matrix metalloproteinase-2-like isoform X1 [Acropora millepora] | XP_029201931.2                                               | This project                                                 |
|                            |                  | MYS              | evm.model.Chr6.535     | myosin heavy chain, striated muscle                          | v1g93924                                                     | https://doi.org/10.1016/j.cell.2018.05.019                   |
|                            |                  | TMP              | evm.model.Chr10.546    | tropomyosin                                                  | /                                                            | https://doi.org/10.1038%2Fnature11180                        |
|                            |                  | MYL12B           | evm.model.Chr4.1556    | myosin regulatory light chain 12B                            | /                                                            | https://doi.org/10.1007%2Fs10162-021-00796-1                 |
|                            |                  | MLC5             | evm.model.Chr9.847     | myosin-2 essential light chain                               | /                                                            | https://doi.org/10.1242/dev.039412                           |
| 11. Progenitor cells       | C14, C18         | NUPR3            | evm.model.Chr10.1218   | nuclear polyadenylated RNA-binding protein 3-like [Acropora  millepora] | XP_029211799.2                                               | This project                                                 |
|                            |                  | P5673            | evm.model.Chr5.88      | hypothetical protein P5673_030106 [Acropora cervicornis]     | KAK2549431.1                                                 | This project                                                 |
|                            |                  | CTHRC1           | evm.model.Chr6.635     | collagen triple helix repeat-containing protein 1-like isoform  X3 [Acropora millepora] | XP_029193101.2                                               | This project                                                 |
|                            |                  | LOC114977096     | evm.model.Chr1.667.2   | uncharacterized protein LOC114977096 isoform X1 [Acropora millepora] | XP_029213539.2                                               | This project                                                 |
|                            |                  | HRJD9            | evm.model.Chr4.73      | JmjC domain containing protein 9                             | v1g172835                                                    | https://doi.org/10.1093/gbe/evz021                           |
|                            |                  | HRJD7            | evm.model.Chr8.1430    | JmjC domain containing protein 7                             | v1g244906                                                    | https://doi.org/10.1093/gbe/evz021                           |
|                            |                  | HRJD3            | evm.model.Chr10.551    | JmjC domain containing protein 3                             | v1g245101                                                    | https://doi.org/10.1093/gbe/evz021                           |
| 12. Immune cells           | C1, C15, C21     | LACE1            | evm.model.Chr8.668     | PREDICTED: lactation elevated protein 1-like [Acropora digitifera] | XP_015772185.1                                               | This project                                                 |
|                            |                  | LOC114961207     | evm.model.Chr7.585     | uncharacterized protein LOC114961207 [Acropora millepora]    | XP_029195664.1                                               | This project                                                 |
|                            |                  | LOC114958720     | evm.model.Chr12.518    | uncharacterized protein LOC114958720 isoform X1 [Acropora millepora] | XP_044178762.1                                               | This project                                                 |
|                            |                  | AIF1L            | evm.model.Chr10.1356   | allograft inflammatory factor 1-like                         | XP_015755194.1                                               | https://doi.org/10.1016/j.fsi.2017.05.063                    |
|                            |                  | IRF1L            | evm.model.Chr7.581     | interferon regulatory factor 1-like                          | XP_022802946.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | IRF2             | evm.model.Chr7.297     | interferon regulatory factor 2                               | XP_022806069.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | IRF2L_1          | evm.model.Chr9.1586    | interferon regulatory factor 2-like                          | XP_022791712.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |
|                            |                  | CTSL1_2          | evm.model.Chr9.1137    | cathepsin L1-like                                            | /                                                            | https://doi.org/10.1016/j.dci.2018.08.014                    |
|                            |                  | IRF2L_2          | evm.model.Chr9.1619    | interferon regulatory factor 2-like                          | XP_022791712.1                                               | https://doi.org/10.1016/j.cell.2021.04.005                   |

## Installation

This script and its dependencies are specifically developed for environments where R scripts are executed, including Linux, Windows, and Mac operating systems, to ensure smooth operation across different platforms. For your convenience, we have provided the following commands, which you can choose from based on your operating system, to download the script.

```shell
git clone https://github.com/coral-seu/coral.git
cd coral/dev_scripts
```

### Example Code

```shell
Rscript seurat_umap_tsne_coral.R
--path AmBD_2vsAmBD_3vsAmBD_4_tSNE.csv
--cluster AmBD_2vsAmBD_3vsAmBD_4_cluster.csv
--name file
--projection tsne
--sample coral_args
--outdir ./
```

### R Requirements

Before running this script, you must ensure that your R environment is version 4.1.3. Additionally, in order to successfully execute the R script, you need to install the following listed dependent packages. These dependencies contain the necessary functions and features required for the script to run.

```R
argparse
ggplot2
tidyverse
pals
ggtext
```

### Input File Requirements Clarification

1. The `path` file should contain the output results of the tSNE  or UMAP algorithm. These results typically include the coordinate information of each cell in the reduced-dimensionality space, facilitating subsequent visualization analysis.
2. The `cluster` file needs to include barcode (unique identifiers for individual cells) and their corresponding cluster information. The barcode is used to uniquely identify a single cell, while the cluster information reflects the categorization of these cells after cluster analysis.
3. The `name` file should provide an extended or explanatory name for each cluster in the `cluster` file. These names help researchers intuitively understand the cell type or state represented by each cluster.

### Detailed explanation of other parameter requirements

1. The `projection` parameter should correspond to whether the input result is from tSNE or UMAP. If the result is from tSNE, the parameter should be set to 'tsne'; if it's from UMAP, it should be set to 'umap'.
2. The `sample` parameter refers to the file name information for the output. This allows users to specify a unique identifier or name for the output files, ensuring organized and recognizable results.
3. The `outdir` parameter corresponds to the directory path where the output files will be saved. This enables users to define a specific folder where all the generated output will be stored, facilitating easy access and management of the results.

### Example output results
![coral_args_Cluster_tSNE_1](https://github.com/coral-seu/coral/assets/172483957/cf0a934e-431c-4886-ab4f-c2dd22df4de8)
![coral_args_Cluster_tSNE_2](https://github.com/coral-seu/coral/assets/172483957/f00b0885-cbba-4417-8837-7f338c2c6c31)
![coral_args_Cluster_tSNE_3](https://github.com/coral-seu/coral/assets/172483957/6a631a36-f006-4272-a785-40d7dd99fba6)
![coral_args_Cluster_tSNE_4](https://github.com/coral-seu/coral/assets/172483957/16b6cc8c-2561-4939-80a4-4b4319df975a)

There are four different types of output results available for users to choose from, including color text, black-and-white text, text with no border, or text annotation disabled for clustering results.
