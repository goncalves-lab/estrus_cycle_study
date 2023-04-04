#Cell-to-cell communication analysis - Visium

#------Libraries--------

library(Seurat) 
library(hash)
library(readr)



#------Functions---------

celltype_signal_map = function(expr_table, ligand_map, receptor_map, prop) {
  #computes the amount of ligand that is captured by the cell type, relative to their expression of the receptor
  signal_map = data.frame(row.names = rownames(expr_table))
  for (i in 1:nrow(mouse_db)) {
    ligand = mouse_db$ligand[i]
    receptor1 = mouse_db$receptor1[i]
    receptor2 = mouse_db$receptor2[i]
    
    if (is.na(receptor2) & ligand%in%colnames(expr_table) & receptor1%in%colnames(expr_table)) {
      receptor = expr_table[,receptor1]
      tot_receptor = receptor_map[rownames(expr_table),receptor1]
      lig = ligand_map[rownames(expr_table), ligand]
      
      signal_map[,paste(ligand, receptor1, sep="_")] = lig * receptor * prop / tot_receptor
    }
    
    if (!is.na(receptor2) & ligand%in%colnames(expr_table) & receptor1%in%colnames(expr_table) & receptor2%in%colnames(expr_table)) {
      receptor = apply(expr_table[,colnames(expr_table)%in%c(receptor1, receptor2)], 1, min, na.rm = T)
      tot_receptor = apply(receptor_map[rownames(expr_table),colnames(receptor_map)%in%c(receptor1, receptor2)], 1, min, na.rm=T)
      lig = ligand_map[rownames(expr_table), ligand]
      
      signal_map[,paste(ligand, receptor1, receptor2, sep="_")] = lig * receptor * prop / tot_receptor
    }
  }
  
  return(signal_map)
}



#------Data----------

load(file="/omics/odcf/analysis/OE0538_projects/DO-0007/mouse_db.Rdata")
pairs_rec = rownames(mouse_db)
pairs_rec = strsplit(pairs_rec, "_")
mouse_db$receptor1 = stringr::str_to_title(sapply(pairs_rec, function(X){return(X[2])}))
mouse_db$receptor2 = stringr::str_to_title(sapply(pairs_rec, function(X){return(X[3])}))
#mouse_db: data.frame with at least these three columns
#                   ligand receptor1 receptor2
# IL1B_IL1R2          Il1b     Il1r2      <NA>
# IL1B_IL1R1_IL1RAP   Il1b     Il1r1    Il1rap
# IL6_IL6RA            Il6     Il6ra      <NA>
# IL6_IL6ST            Il6     Il6st      <NA>


#------Main---------

cell_prop_hash = hash()
APC_expr_hash = hash()
BC_expr_hash = hash()
EC_expr_hash = hash()
GC_expr_hash = hash()
LC_expr_hash = hash()
SC_expr_hash = hash()
TC_expr_hash = hash()
EpC_expr_hash = hash()

total_expr_hash = hash()
ligand_map_hash = hash()
receptor_map_hash = hash()
signal_SC_map_hash = hash()
signal_APC_map_hash = hash()
signal_BC_map_hash = hash()
signal_EC_map_hash = hash()
signal_EpC_map_hash = hash()
signal_GC_map_hash = hash()
signal_LC_map_hash = hash()
signal_TC_map_hash = hash()





for (slide in c("V11M25-311_B1", "V11M25-312_A1", "V11M25-312_B1", "V11M25-312_C1", "V11M25-312_D1",  #ovary 3M
                "V11M25-373_C1", "V11M25-373_D1", "V19D02-108_A1", "V19D02-108_B1", "V19D02-108_C1", "V19D02-108_D1" #ovary 18M
)) {
  
  if (! slide%in%names(signal_BC_map_hash)) {
    setwd("/dir/where/the/data/is/")
    coords = read.csv(paste0(slide, "tissue_positions_list.csv"), header=F, row.names = 1) #tissue position list in the same format as outputed by spaceranger
    colnames(coords) = c("tissue", "row", "col", "imagerow", "imagecol")
    coords$imagerow = -coords$imagerow
    
    cell_prop = read.csv(paste0(slide, "_corrected_cell_proportions.csv"), row.names = 1) #cell type proportions from DestVI deconvolution
    
    #Cell-type specific gene expression profiles outputed by DestVI deconvolution
    APC_expr = read.table(paste0(slide, "_APC_gene_expression_thresholded.csv"), row.names = 1, header = T, sep = "\t")
    BC_expr = read.table(paste0(slide, "_BC_gene_expression_thresholded.csv"), row.names = 1, header = T, sep = "\t")
    EC_expr = read.table(paste0(slide, "_EC_gene_expression_thresholded.csv"), row.names = 1, header = T, sep = "\t")
    GC_expr = read.table(paste0(slide, "_GC_gene_expression_thresholded.csv"), row.names = 1, header = T, sep = "\t")
    LC_expr = read.table(paste0(slide, "_LC_gene_expression_thresholded.csv"), row.names = 1, header = T, sep = "\t")
    SC_expr = read.table(paste0(slide, "_SC_gene_expression_thresholded.csv"), row.names = 1, header = T, sep = "\t")
    TC_expr = read.table(paste0(slide, "_TC_gene_expression_thresholded.csv"), row.names = 1, header = T, sep = "\t")
    EpC_expr = read.table(paste0(slide, "_EpC_gene_expression_thresholded.csv"), row.names = 1, header = T, sep = "\t")
    
    #making sure that all spots are properly deconvoluted
    cell_prop = cell_prop[rowSums(cell_prop, na.rm = T)>0.99, ]
    
    #distance matrix between spots, and selecting neighborhoods spots
    dist_mat = as.matrix(dist(coords[rownames(cell_prop),c("imagerow", "imagecol")]))
    neighborhood = dist_mat < min(dist_mat[dist_mat!=0])*1.1

    #computing total expression of each spot by summing the contributions of each cell type
    total_expr = as.data.frame(matrix(data = 0, nrow = nrow(cell_prop), ncol = ncol(APC_expr)))
    colnames(total_expr) = colnames(APC_expr)
    rownames(total_expr) = rownames(cell_prop)
    total_expr[rownames(APC_expr),] = total_expr[rownames(APC_expr),] + APC_expr * cell_prop[rownames(APC_expr), "APC"]
    total_expr[rownames(BC_expr),] = total_expr[rownames(BC_expr),] + BC_expr * cell_prop[rownames(BC_expr), "BC"]
    total_expr[rownames(EC_expr),] = total_expr[rownames(EC_expr),] + EC_expr * cell_prop[rownames(EC_expr), "EC"]
    total_expr[rownames(EpC_expr),] = total_expr[rownames(EpC_expr),] + EpC_expr * cell_prop[rownames(EpC_expr), "EpC"]
    total_expr[rownames(GC_expr),] = total_expr[rownames(GC_expr),] + GC_expr * cell_prop[rownames(GC_expr), "GC"]
    total_expr[rownames(LC_expr),] = total_expr[rownames(LC_expr),] + LC_expr * cell_prop[rownames(LC_expr), "LC"]
    total_expr[rownames(SC_expr),] = total_expr[rownames(SC_expr),] + SC_expr * cell_prop[rownames(SC_expr), "SC"]
    total_expr[rownames(TC_expr),] = total_expr[rownames(TC_expr),] + TC_expr * cell_prop[rownames(TC_expr), "TC"]
    
    
    #generating the ligand availability map by spreading it over the neighborhood
    ligand_map = data.frame(row.names = rownames(cell_prop))
    
    for (ligand in unique(mouse_db$ligand)) {
      vect = rep(0, nrow(neighborhood))
      names(vect) = rownames(neighborhood)
      
      if (ligand%in%colnames(total_expr)) {
        for (spot in rownames(neighborhood)) {
          neighbs = colnames(neighborhood)[neighborhood[spot,]==T]
          vect[neighbs] = vect[neighbs] + total_expr[spot,ligand]/length(neighbs)
        }
        
        ligand_map[,ligand] = vect
      }
    }
    
    
    #generating the receptor map of how much receptor is expressed (no matter the cell type)
    receptor_map = data.frame(row.names = rownames(total_expr))
    for (i in 1:nrow(mouse_db)) {
      ligand = mouse_db$ligand[i]
      receptor1 = mouse_db$receptor1[i]
      receptor2 = mouse_db$receptor2[i]
      
      if (is.na(receptor2) & receptor1%in%colnames(total_expr)) {
        tot_receptor = total_expr[,receptor1]
        
        receptor_map[,paste(ligand, receptor1, sep="_")] = tot_receptor
      }
      
      if (!is.na(receptor2) & ligand%in%colnames(total_expr) & receptor1%in%colnames(total_expr) & receptor2%in%colnames(total_expr)) {
        tot_receptor = apply(total_expr[,colnames(total_expr)%in%c(receptor1, receptor2)], 1, min, na.rm=T)
        
        receptor_map[,paste(ligand, receptor1, receptor2, sep="_")] = tot_receptor
      }
    }
    
    
    #computing amount of ligands captured by each cell_type
    signal_SC_map = celltype_signal_map(SC_expr, ligand_map, total_expr, cell_prop[rownames(SC_expr),"SC"])
    signal_APC_map = celltype_signal_map(APC_expr, ligand_map, total_expr, cell_prop[rownames(APC_expr),"APC"])
    signal_BC_map = celltype_signal_map(BC_expr, ligand_map, total_expr, cell_prop[rownames(BC_expr),"BC"])
    signal_EC_map = celltype_signal_map(EC_expr, ligand_map, total_expr, cell_prop[rownames(EC_expr),"EC"])
    signal_EpC_map = celltype_signal_map(EpC_expr, ligand_map, total_expr, cell_prop[rownames(EpC_expr),"EpC"])
    signal_GC_map = celltype_signal_map(GC_expr, ligand_map, total_expr, cell_prop[rownames(GC_expr),"GC"])
    signal_LC_map = celltype_signal_map(LC_expr, ligand_map, total_expr, cell_prop[rownames(LC_expr),"LC"])
    signal_TC_map = celltype_signal_map(TC_expr, ligand_map, total_expr, cell_prop[rownames(TC_expr),"TC"])
    
    
    #put the results in the hash
    cell_prop_hash[slide] = cell_prop
    APC_expr_hash[slide] = APC_expr
    BC_expr_hash[slide] = BC_expr
    EC_expr_hash[slide] = EC_expr
    GC_expr_hash[slide] = GC_expr
    LC_expr_hash[slide] = LC_expr
    SC_expr_hash[slide] = SC_expr
    TC_expr_hash[slide] = TC_expr
    EpC_expr_hash[slide] = EpC_expr
    total_expr_hash[slide] = total_expr
    ligand_map_hash[slide] = ligand_map
    receptor_map_hash[slide] = receptor_map
    signal_SC_map_hash[slide] = signal_SC_map
    signal_APC_map_hash[slide] = signal_APC_map
    signal_BC_map_hash[slide] = signal_BC_map
    signal_EC_map_hash[slide] = signal_EC_map
    signal_EpC_map_hash[slide] = signal_EpC_map
    signal_GC_map_hash[slide] = signal_GC_map
    signal_LC_map_hash[slide] = signal_LC_map
    signal_TC_map_hash[slide] = signal_TC_map
    
  }
}




#Computing ligand-receptor interaction scores for stromal cells based on the amount of ligand they capture and their expression of the receptor
receptor_map_SC_hash = hash()
for (slide in c("V11M25-311_B1", "V11M25-312_A1", "V11M25-312_B1", "V11M25-312_C1", "V11M25-312_D1"
                ,"V11M25-373_C1", "V11M25-373_D1", "V19D02-108_A1", "V19D02-108_B1", "V19D02-108_C1", "V19D02-108_D1"
)) {
  receptor_map = data.frame(row.names = rownames(SC_expr_hash[[slide]]))
  for (i in 1:nrow(mouse_db)) {
    ligand = mouse_db$ligand[i]
    receptor1 = mouse_db$receptor1[i]
    receptor2 = mouse_db$receptor2[i]
    
    if (is.na(receptor2) & receptor1%in%colnames(SC_expr_hash[[slide]])) {
      tot_receptor = SC_expr_hash[[slide]][,receptor1]
      
      receptor_map[,paste(ligand, receptor1, sep="_")] = tot_receptor
    }
    
    if (!is.na(receptor2) & ligand%in%colnames(SC_expr_hash[[slide]]) & receptor1%in%colnames(SC_expr_hash[[slide]]) & receptor2%in%colnames(SC_expr_hash[[slide]])) {
      tot_receptor = apply(SC_expr_hash[[slide]][,colnames(SC_expr_hash[[slide]])%in%c(receptor1, receptor2)], 1, min, na.rm=T)
      
      receptor_map[,paste(ligand, receptor1, receptor2, sep="_")] = tot_receptor
    }
  }
  receptor_map_SC_hash[slide] = receptor_map
}

all_LR_slides = data.frame(Row.names = colnames(signal_SC_map_hash[["V11M25-311_B1"]]), V1 = colMeans(signal_SC_map_hash[["V11M25-311_B1"]]*
                                                                                                        receptor_map_SC_hash[["V11M25-311_B1"]][rownames(signal_SC_map_hash[["V11M25-311_B1"]]), colnames(signal_SC_map_hash[["V11M25-311_B1"]])]))
for (slide in c("V11M25-312_A1", "V11M25-312_B1", "V11M25-312_C1", "V11M25-312_D1"
                ,"V11M25-373_C1", "V11M25-373_D1", "V19D02-108_A1", "V19D02-108_B1", "V19D02-108_C1", "V19D02-108_D1"
)) {
  all_LR_slides = merge(x = all_LR_slides,
                        y = data.frame(Row.names = colnames(signal_SC_map_hash[[slide]]), 
                                       V1 = colMeans(signal_SC_map_hash[[slide]]*
                                                       receptor_map_SC_hash[[slide]][rownames(signal_SC_map_hash[[slide]]), colnames(signal_SC_map_hash[[slide]])]
                                       )),
                        by = "Row.names", all = T)
}
rownames(all_LR_slides) = all_LR_slides$Row.names
all_LR_slides = all_LR_slides[,-1]
colnames(all_LR_slides) = c("V11M25-311_B1", "V11M25-312_A1", "V11M25-312_B1", "V11M25-312_C1", "V11M25-312_D1",
                            "V11M25-373_C1", "V11M25-373_D1", "V19D02-108_A1", "V19D02-108_B1", "V19D02-108_C1", "V19D02-108_D1")

