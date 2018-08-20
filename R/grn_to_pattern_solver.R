options(warn=-1)

################
## FORMAT     ##
################

#' Format the resulting GRN 
#' as shown in the observations.spec file
#'
#' @param vgene_ls a pair of gene name, expression level output
#' @param formula character string that describes the logical formula
#'        which output is vgene_ls[2]
formatting <- function(vgene_ls, formula) {
	res <- paste0("\t$", vgene_ls[1], " := {\n", 
			"\t\t$formula := {\n", 
			"\t\t\t", formula, ";\n",
		 	"\t\t\t$output := ", vgene_ls[2], ";\n", 
			"\t\t};\n",
			"\t}")
	return(res)
}

################
## TOOLS      ##
################

#' @param gene original string: gene = "GENE-1"
#' @return vgene_ls where vgene_ls[1] = "Gene" and vgene_ls[2] = gene value = 1
get_vgene_value <- function(gene) {
	tmp = strsplit(gene, "-")[[1]]
	vgene = tolower(paste0(tmp[1:(length(tmp)-1)], collapse="-"))
	vgene = paste0(toupper(substr(vgene, 1, 1)), substr(vgene, 2, nchar(vgene)))
	value = as.integer(tmp[length(tmp)])
	if (is.na(value)) value <- 1
	return(c(vgene, value))
}

#' @param element character string like "[~][gene name]-[binary value]"
#' @return res character string equivalent to input, like "[gene name]=[integer value]"
#'        or "( [gene name]=[integer value] or [gene name]=[integer value] or ... )"
to_multi <- function(element) {
	## Ignore KO-related variables
	if (grepl("^~KO[_]", element)) return(NULL)
	else {
		## Ignore FE-related variables
		if (grepl("^FE[_]", element)) return(NULL)
		else {
			## Test if expression is a negation
			negative <- grepl("^~", element)
			## Get gene name + value
			if (negative) gene <- strsplit(element, "~")[[1]][2]
			else gene <- element
			vgene_ls <- get_vgene_value(gene)
			## If negation, then allow all lower values
			## i.e. if element = "~Gene-i"
			## then res = "( Gene=0 or Gene=1 or ... or Gene=i-1 )"
			if (negative) {
				return(paste0("( ", 
					paste0(sapply(0:(as.integer(vgene_ls[2])-1), function(v) {
						paste0(vgene_ls[1], "=", v)
					}), collapse=" or "), " )"))
			}
			## Otherwise, if element="Gene-i", res = "Gene=i"
			else return(paste0(vgene_ls[1], "=", vgene_ls[2]))
		}
	}
}

################
## MAIN       ##
################

#' Turn the GRF result file (from the previous GRN solver)
#' into a readable string for the PATTERN solver
#'
#' @param filename character string of the (relative path) + file name of the
#'                GRF result file (without the .txt extension)
#' @return None (writes the resulting string into a new file)
grn_to_pattern_solver <- function(filename) {
	new_filename = paste0(filename, "_converted.txt")
	con = file(paste0(filename, ".txt"), "r")
	lines = readLines(con)
	close(con)
	str = unlist(sapply(lines, function(line) {
		## Ignore KO-related variables
		if (grepl("GRF[(]KO[_]", line)) {NULL}
		else {
			## Ignore FE-related variables
			if (grepl("GRF[(]FE[_]", line)) {NULL}
			else {
				## Get GRF line
				tmp = strsplit(strsplit(line, "GRF[(]")[[1]][2], "[)] = ")[[1]]
				## Get target gene + output value of formula
				vgene_ls = get_vgene_value(tmp[1])
				## Get raw formula and split it into disjunctive clauses
				## then into conjunctive clauses (of one litteral)
				formula = lapply(strsplit(tmp[2], " [+] ")[[1]], 
					function(e) strsplit(e, "[*]")[[1]])
				## Convert litterals into equalities (with integer values)
				## and rebuild conjunctive clauses
				formula = unlist(sapply(formula, function(clause) {
					values <- unlist(sapply(clause, to_multi))
					values <- values[!is.null(values)]
					if (!is.null(values)) {
						paste0("( ", paste0(values, collapse=" and "), " )")
					}
				}))
				formula = formula[formula != "( NULL )"]
				## Rebuild disjunctive clauses
				formula = paste0("( ", paste0(formula, collapse=" or "), " )")
				## Final formatting of the GRF
				formatting(vgene_ls, formula)
			}
		}
	}))
	str = paste0(str[!is.null(str)], collapse="\n")
	## Write resulting character string to another file
	write(str, file=new_filename)
	print("MSG: Conversion done!")
}
