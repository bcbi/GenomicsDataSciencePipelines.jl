using DataFrames
using CSV

data_dir =
prior_file = 
variant_file =

work_dir =

function sorted_counts(df, col)
	sort(combine(groupby(df, col), nrow), :nrow, rev=true)
end

vdf = CSV.read(joinpath(data_dir, variant_file), DataFrame)
CSV.write(joinpath(work_dir, "accessions.csv"), accessions)
accessions = sorted_counts
filter(x->x.nrow > 1, accessions)
vdf.IDL_City |> unique |> sort .|> println;

