#' Predefined TLS gene signatures
#'
#' @description
#' A list of example human TLS gene signatures included with the package. Users can supply
#' their own signatures to the scoring functions.
#'
#' 12-Gene Signature Source- Coppola D, Nebozhyn M, Khalil F, Dai H, Yeatman T, Loboda A, Mulé JJ. 
#' Unique ectopic lymph node-like structures present in human primary colorectal carcinoma are identified by immune gene array profiling. 
#' Am J Pathol. 2011 Jul;179(1):37-45. doi: 10.1016/j.ajpath.2011.03.007. Epub 2011 May 3. PMID: 21703392; PMCID: PMC3123872.
#'
#' TLS Imprint Source- Meylan M, Petitprez F, Becht E, Bougoüin A, Pupier G, Calvez A, Giglioli I, Verkarre V, Lacroix G, Verneau J, Sun CM, Laurent-Puig P, Vano YA, Elaïdi R, Méjean A, Sanchez-Salas R, Barret E, Cathelineau X, Oudard S, Reynaud CA, de Reyniès A, Sautès-Fridman C, Fridman WH. 
#' Tertiary lymphoid structures generate and propagate anti-tumor antibody-producing plasma cells in renal cell cancer. 
#' Immunity. 2022 Mar 8;55(3):527-541.e5. doi: 10.1016/j.immuni.2022.02.001. Epub 2022 Feb 28. PMID: 35231421.
#'
#' 9-Gene Signature Source- Cabrita R, Lauss M, Sanna A, Donia M, Skaarup Larsen M, Mitra S, Johansson I, Phung B, Harbst K, Vallon-Christersson J, van Schoiack A, Lövgren K, Warren S, Jirström K, Olsson H, Pietras K, Ingvar C, Isaksson K, Schadendorf D, Schmidt H, Bastholt L, Carneiro A, Wargo JA, Svane IM, Jönsson G. 
#' Tertiary lymphoid structures improve immunotherapy and survival in melanoma. 
#' Nature. 2020 Jan;577(7791):561-565. doi: 10.1038/s41586-019-1914-8. Epub 2020 Jan 15. Erratum in: Nature. 2020 Apr;580(7801):E1. doi: 10.1038/s41586-020-2155-6. PMID: 31942071.
#'
#' 36-Gene Signature Source- Mughal SS, Reiss Y, Felsberg J, Meyer L, Macas J, Schlue S, Starzetz T, Köhrer K, Fehm T, Müller V, Lamszus K, Schadendorf D, Helfrich I, Wikman H, Berghoff A, Brors B, Plate KH, Reifenberger G. 
#' Identification and characterization of tertiary lymphoid structures in brain metastases. 
#' Acta Neuropathol Commun. 2025 May 3;13(1):91. doi: 10.1186/s40478-025-02007-x. PMID: 40319321; PMCID: PMC12049775.
#'
#' @format A named list of character vectors; each element is a vector of gene symbols.
#' @examples
#' human_tls_signatures$TLS_12_GENE_SIGNATURE
#' @export
human_tls_signatures <- list(
  TLS_12_GENE_SIGNATURE = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",  "CXCL9", "CXCL10", "CXCL11", "CXCL13"),
  TLS_IMPRINT = c("IGHA1", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM", "IGKC", "IGLC1", "IGLC2", "IGLC3", "JCHAIN", "CD79A", "FCRL5", "MZB1", "SSR4", "XBP1", "TRBC2", "Il7R", "CXCL12", "LUM", "C1QA", "C7", "CD52", "APOE", "PLTP", "PTGDS", "PIM2", "DERL3"),
  TLS_9_GENE_SIGNATURE = c("CD79B", "CD1D", "CCR6", "LAT", "SKAP1", "CETP", "EIF1AY", "RBP5", "PTGDS"),
  TLS_36_GENE_SIGNATURE = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CD79B", "CETP", "LAT", "CD1D", "PTGDS", "CXCR5", "SELL", "ICOS", "SH2D1A", "TIGIT", "PDCD1", "BANK1", "CD22", "CD79A", "CR1", "FCRl2", "MS4A1", "FCER2", "LAMP3", "CD86", "CD80", "CD83", "CCR7")
)