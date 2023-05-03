using ConfParser
using CSV
using DataFrames
using Dates
using UnicodePlots

# Unique values in df.col, sorted by number of occurrences
function sorted_counts(df, col)
	sort(combine(groupby(df, col), nrow), :nrow, rev=true)
end

"""
# Find data dictionary entry for a column
# data_dictionary("Conditions") # Example
"""
function data_dictionary(x)
	ddf[ddf[!,"Variable Name"] .== uppercase(x),:]
end

# Read in options
conf = ConfParse("./config.ini")
parse_conf!(conf)
user_vars = [ "data_dir", "prior_file", "variant_file", "work_dir", "dict_file", "accession_file"]
for var in user_vars
	Meta.parse("""$var = retrieve(conf, "local", "$var")""") |> eval
end

# Read in data
vdf = CSV.File(joinpath(data_dir, variant_file),
               dateformat="mm/dd/yyyy",
               types=Dict(:idl_zipcode=>String, test_site_zip=>String) ) |> DataFrame
ddf = CSV.read(joinpath(work_dir, dict_file), DataFrame, select=1:6)
gdf = CSV.read(joinpath(work_dir, accession_file), DataFrame)

# Remove unneeded columns
#select!(ddf, Not(:Column7)) # or just read in the first 6 columns only
# Columns where all values are the same
bad_columns = names(vdf)[ vdf |> eachcol .|> allequal |> findall ]
# Columns where all values are missing
#names(vdf)[ vdf |> eachcol .|> skipmissing .|> isempty |> findall ]
# Delete these columns
select!(vdf, Not( bad_columns ))

# Set proper datatype for true/false survey data
for col in names(vdf)
	all_values = unique(vdf[!,col])
	if all_values == [0,1]
		vdf[!,col] = Bool.(vdf[!,col])
	elseif all_values == ["No", "Yes"]
		d = Dict(["No", "Yes"] .=> [false,true])
		vdf[!,col] = vdf[!,col] .|> x->d[x]
	end
end

# Fix names
rename!(vdf, names(vdf) .=> uppercase.(names(vdf)))
ddf[!,"Variable Name"] = uppercase.(ddf[!,"Variable Name"])

# Summaries
#show(describe(vdf), allrows=true, allcols=true)

# Ad-hoc checks
#accessions = sorted_counts(
#CSV.write(joinpath(work_dir, "accessions.csv"), accessions)
#filter(x->x.nrow > 1, accessions)
#vdf.IDL_City |> unique |> sort .|> println;
#vdf.IDL_age |> histogram
#vdf.IDL_age |> boxplot
#filter(x->occursin("0",x), vdf.Accession_Number |> skipmissing) |> unique |> length
#vdf.IDL_specimen_collection_date |> skipmissing |> extrema
#x = names(vdf);
#y = ddf[!, "Variable Name"];
#x ∩ y # Columns found in the data dictionary
#x[x.∉Ref(y)] # Columns not found in the data dictionary
#y[y.∉Ref(x)] # Dictionary entries not found in the variants columns

# Connect accession numbers with data set
m = match.(r"^hCoV-19/USA/(.*)-(.*)-(.*)/([0-9]{4})$", gdf.strain)
n = unique(filter(x->!isnothing(x),m) .|> x->x.captures[3])
n ∩ vdf.ACCESSION_NUMBER
