
# Equivalent competition
nci_eq <- function(neighbors, alpha, beta){
  raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
  return(sum(raw))
}

# Intraspecific vs. interspecific competition
nci_int <- function(neighbors, alpha, beta, intra, inter){
  raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
  cons <- which(neighbors$sps_comp == focal_sps)
  hets <- which(neighbors$sps_comp != focal_sps)
  nci_con <- raw[cons] * intra
  nci_het <- raw[hets] * inter
  return(sum(c(nci_con, nci_het)))
}

# Species-specific competition (4 competitor species)
nci_ss4 <- function(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4){
  raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
  nci1 <- raw[which(neighbors$sps_comp == comps[1])] * lmd1
  nci2 <- raw[which(neighbors$sps_comp == comps[2])] * lmd2
  nci3 <- raw[which(neighbors$sps_comp == comps[3])] * lmd3
  nci4 <- raw[which(neighbors$sps_comp == comps[4])] * lmd4
  return(sum(c(nci1, nci2, nci3, nci4)))
}

# Species-specific competition (5 competitor species)
nci_ss5 <- function(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4, lmd5){
  raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
  nci1 <- raw[which(neighbors$sps_comp == comps[1])] * lmd1
  nci2 <- raw[which(neighbors$sps_comp == comps[2])] * lmd2
  nci3 <- raw[which(neighbors$sps_comp == comps[3])] * lmd3
  nci4 <- raw[which(neighbors$sps_comp == comps[4])] * lmd4
  nci5 <- raw[which(neighbors$sps_comp == comps[5])] * lmd5
  return(sum(c(nci1, nci2, nci3, nci4, nci5)))
}

# Species-specific competition (8 competitor species)
nci_ss8 <- function(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4, lmd5,
                lmd6, lmd7, lmd8){
  raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
  nci1 <- raw[which(neighbors$sps_comp == comps[1])] * lmd1
  nci2 <- raw[which(neighbors$sps_comp == comps[2])] * lmd2
  nci3 <- raw[which(neighbors$sps_comp == comps[3])] * lmd3
  nci4 <- raw[which(neighbors$sps_comp == comps[4])] * lmd4
  nci5 <- raw[which(neighbors$sps_comp == comps[5])] * lmd5
  nci6 <- raw[which(neighbors$sps_comp == comps[6])] * lmd6
  nci7 <- raw[which(neighbors$sps_comp == comps[7])] * lmd7
  nci8 <- raw[which(neighbors$sps_comp == comps[8])] * lmd8
  return(sum(c(nci1, nci2, nci3, nci4, nci5, nci6, nci7, nci8)))
}

# Species-specific competition (9 competitor species)
nci_ss9 <- function(neighbors, alpha, beta, lmd1, lmd2, lmd3, lmd4, lmd5,
                lmd6, lmd7, lmd8, lmd9){
  raw <- (neighbors$dbh_comp ^ alpha) / (neighbors$prox ^ beta)
  nci1 <- raw[which(neighbors$sps_comp == comps[1])] * lmd1
  nci2 <- raw[which(neighbors$sps_comp == comps[2])] * lmd2
  nci3 <- raw[which(neighbors$sps_comp == comps[3])] * lmd3
  nci4 <- raw[which(neighbors$sps_comp == comps[4])] * lmd4
  nci5 <- raw[which(neighbors$sps_comp == comps[5])] * lmd5
  nci6 <- raw[which(neighbors$sps_comp == comps[6])] * lmd6
  nci7 <- raw[which(neighbors$sps_comp == comps[7])] * lmd7
  nci8 <- raw[which(neighbors$sps_comp == comps[8])] * lmd8
  nci9 <- raw[which(neighbors$sps_comp == comps[9])] * lmd9
  return(sum(c(nci1, nci2, nci3, nci4, nci5, nci6, nci7, nci8, nci9)))
}