using CSV
using DataFrames
using Dates
using UnicodePlots

data_dir =
prior_file = 
variant_file =
work_dir =
dict_file =
accession_file = 


function sorted_counts(df, col)
	sort(combine(groupby(df, col), nrow), :nrow, rev=true)
end

# Find data dictionary entry for a column
function data_dictionary(x)
	ddf[ddf[!,"Variable Name"] .== uppercase(x),:]
end

vdf = CSV.File(joinpath(data_dir, variant_file), dateformat="mm/dd/yyyy") |> DataFrame
CSV.write(joinpath(work_dir, "accessions.csv"), accessions)
accessions = sorted_counts
filter(x->x.nrow > 1, accessions)
vdf.IDL_City |> unique |> sort .|> println;
vdf.IDL_age |> histogram
vdf.IDL_age |> boxplot

filter(x->occursin("0",x), vdf.Accession_Number |> skipmissing) |> unique |> length

Dates.DateTime.(skipmissing(vdf.IDL_specimen_collection_date), dateformat"mm/dd/yyyy") |> extrema

ddf = CSV.read(joinpath(work_dir, dict_file), DataFrame)
select!(ddf, Not(:Column7)) # or just read in the first 6 columns only

x = names(vdf);
y = ddf[!, "Variable Name"];

x ∩ y # Columns found in the data dictionary
x[x.∉Ref(y)] # Columns not found in the data dictionary
y[y.∉Ref(x)] # Dictionary entries not found in the variants columns

# Fix names
rename!(vdf, names(vdf) .=> uppercase.(names(vdf)))
ddf[!,"Variable Name"] = uppercase.(ddf[!,"Variable Name"])

data_dictionary("Conditions") # Example

cdf=DataFrame
for col in names(vdf)
	values = vdf[!,col] |> skipmissing |> extrema
	isempty(values) & 0
end

# Summaries
show(describe(vdf), allrows=true, allcols=true)

# Columns where all values are the same
bad_columns = names(vdf)[ vdf |> eachcol .|> allequal |> findall ]
		
# Columns where all values are missing
names(vdf)[ vdf |> eachcol .|> skipmissing .|> isempty |> findall ]

# Delete these columns
select!(ddf, Not( bad_columns ))

#TODO: Some Int64 columns should be Bool

# Connect accession numbers with data set
gdf = CSV.read(joinpath(work_dir, accession_file), DataFrame)

		
