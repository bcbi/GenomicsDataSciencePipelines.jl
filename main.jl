using DataFrames
using CSV

data_dir =
prior_file = 
variant_file =

work_dir =

vdf = CSV.read(joinpath(data_dir, variant_file), DataFrame, select=[:Accession_Number])
accessions = sort(combine(groupby(vdf, :Accession_Number), nrow), :nrow, rev=true)
CSV.write(joinpath(work_dir, "accessions.csv"), accessions)

