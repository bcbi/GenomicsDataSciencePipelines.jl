using CSV
using DataFrames
using Dates
using UnicodePlots

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
vdf.IDL_age |> histogram
vdf.IDL_age |> boxplot

filter(x->occursin("0",x), vdf.Accession_Number |> skipmissing) |> unique |> length


Dates.DateTime.(skipmissing(vdf.IDL_specimen_collection_date), "mm/dd/yyyy") |> extrema
