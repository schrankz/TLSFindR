#' Predefined TLS gene signatures
#'
#' @description
#' A list of example mouse TLS gene signatures (converted from human orthologs) included with the package. Users can supply
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
#' mouse_tls_signatures$TLS_12_GENE_SIGNATURE
#' @export
mouse_tls_signatures <- list(
  TLS_12_GENE_SIGNATURE = c("Ccl12",  "Ccl2",   "Ccl3",   "Ccl4",   "Ccl5",   "Ccl11",  "Ccl8",   "Ccl19",  "Ccl21d", "Ccl21b", "Ccl21f", "Ccl21e", "Ccl21a", "Cxcl9",  "Cxcl10", "Cxcl11", "Cxcl13"),
  TLS_IMPRINT = c("Igha", "Ighg2c", "Ighg1",  "Ighg3",  "Ighg2b", "Ighg",   "Ighg2a", "Ighm",   "Igkc",   "Iglc1",  "Iglc2",  "Iglc4", "Jchain", "Cd79a",  "Fcrl5",  "Mzb1",   "Ssr4",   "Xbp1",   "Trbc2",  "Trbc1",  "Cxcl12", "Lum",    "C1qa",   "C7",    "Cd52",   "Apoe",   "Pltp",  "Pim2",   "Derl3"),
  TLS_9_GENE_SIGNATURE = c("Cd79b",    "Cd1d2",    "Cd1d1",    "Ccr6",     "Lat",      "Skap1",    "Eif1a",    "Eif1ad6",  "Eif1ad16", "Eif1ad4", "Eif1ad17", "Eif1ad14", "Eif1ad15", "Eif1ad13", "Eif1ad11", "Eif1ad2",  "Eif1ad7",  "Eif1ad18", "Eif1ad8",  "Eif1ad3", "Eif1ad19", "Eif1ad12", "Eif1ax",   "Ptgds"),
  TLS_36_GENE_SIGNATURE = c("Ccl12",  "Ccl2",   "Ccl3",   "Ccl4",   "Ccl5",   "Ccl11",  "Ccl8",   "Ccl19",  "Ccl21d", "Ccl21b", "Ccl21f", "Ccl21e", "Ccl21a", "Cxcl9",  "Cxcl10", "Cxcl11", "Cxcl12", "Cxcl13", "Cd79b",  "Lat",    "Cd1d2",  "Cd1d1",  "Ptgds",  "Cxcr5", "Sell",   "Icos",   "Sh2d1a", "Tigit",  "Pdcd1",  "Bank1",  "Cd22",   "Cd79a",  "Cr1l",   "Ms4a1",  "Fcer2a", "Lamp3", "Cd86",   "Cd80",   "Cd83",   "Ccr7")
)