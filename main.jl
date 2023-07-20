
#===========================================================

The easiest way to push changes to the secure enclave:
pbcopy < ./main.jl # On Mac, copy the file to the clipboard
# On the server, delete the old file contents and paste in the new

USAGE: julia> include("./main.jl")

Everything is in global scope, so you can
work directly with all the data afterwards.

===========================================================#

#-----------------------------------------------------------
@info "Loading packages"
#-----------------------------------------------------------

using ConfParser
using CSV
using DataFrames
using DataFramesMeta
using Dates
using GMT
using StatsBase
using UnicodePlots
using Plots

#-----------------------------------------------------------
@info "Reading in data from files"
#-----------------------------------------------------------

# Read in options
conf = ConfParse("./config.ini")
parse_conf!(conf)
data_dir = retrieve(conf, "local", "data_dir")
prior_file = retrieve(conf, "local", "prior_file")
variant_file = retrieve(conf, "local", "variant_file")

work_dir = retrieve(conf, "local", "work_dir")
dict_file = retrieve(conf, "local", "dict_file")
accession_file = retrieve(conf, "local", "accession_file")

# Read in data
vdf = CSV.File(joinpath(data_dir, variant_file),
	       dateformat="mm/dd/yyyy",
	       types=Dict(:idl_zipcode=>String, :test_site_zip=>String) ) |> DataFrame
ddf = CSV.read(joinpath(work_dir, dict_file), DataFrame, select=1:6)
gdf = CSV.read(joinpath(work_dir, accession_file), DataFrame)

#-----------------------------------------------------------
@info "Fixing data in variant file"
#-----------------------------------------------------------

# Fix ZIP Codes
# This column is a write-in field.
# Most have 4 digits, but one has 3. None begin wtith a zero.
# I am padding all of them with zeros.
# The other ZIP Code field is fine. All present values have 5 digits
for i in 1:nrow(vdf)
	if !ismissing(vdf[i,:idl_zipcode])
		vdf[i,:idl_zipcode] = lpad(vdf[i,:idl_zipcode], 5, '0')
	end

	# I'm also using this loop to fix a typo in one accession number (an extra dash)
	# UPDATE: Nevermind, this doesn't end up in the joined dataset anyway
	#if !ismissing(vdf[i,:Accession_Number])
	#	vdf[i,:Accession_Number] = replace(vdf[i,:Accession_Number], "--" => "-")
	#end
end

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

# Use Int for age (All values are a whole number of years)
vdf.IDL_age = Int.(vdf.IDL_age)

#-----------------------------------------------------------
@info "Making improvements to datasets"
#-----------------------------------------------------------

# Delete columns where all values are the same
for df in [vdf, gdf]
	bad_columns = names(df)[ df |> eachcol .|> allequal |> findall ]
	for col in bad_columns
		println("Deleting column $col; all values are: $(df[1,col])")
	end
	select!(df, Not( bad_columns ))
end

# For consistency between data sets, use all caps for column names
rename!(vdf, names(vdf) .=> uppercase.(names(vdf)))
rename!(gdf, names(gdf) .=> uppercase.(names(gdf)))
ddf[!,"Variable Name"] = uppercase.(ddf[!,"Variable Name"])

#-----------------------------------------------------------
@info "Connecting both datasets"
#-----------------------------------------------------------

# Add ACCESSION_NUMBER to GISAID data matching the COVID-19 accession IDs
gdf.ACCESSION_NUMBER = let
	m = match.(r"^hCoV-19/USA/RI-RISHL-(.*)/20([0-9]{2})$",gdf.STRAIN)
	x = convert(Vector{Union{Nothing,String}}, fill(nothing,nrow(gdf)))
	for i in 1:nrow(gdf)
		a=m[i]
		isnothing(a) && continue
		x[i] = "RI-" * a[2] * "-" * a[1]
	end
	x;
end

# Join datasets on matching ACCESSION_NUMBER
jdf = innerjoin(vdf, gdf, on=:ACCESSION_NUMBER, matchmissing=:notequal)

# One row of gdf seems to have AGE and SEX interchanged.
# However, it doesn't appear in the innerjoin, so it is ignored.

#-----------------------------------------------------------
@info "Preparing data summary"
#-----------------------------------------------------------

# Summary
q = describe(jdf)
q.nunique = map(x->length(jdf[!,x] |> unique), names(jdf))
show(q, allrows=true, allcols=true)
println()

#-----------------------------------------------------------
@info "Drawing Figures"
#-----------------------------------------------------------

Plots.savefig(Plots.histogram(jdf.IDL_SPECIMEN_COLLECTION_DATE, bins=80), "collection_dates.png")

# Geographical plot
cpop = CSV.File("/data/ursa_research/ris3/Ashlin/zips.csv", types=Dict(:ZIP=>String)) |> DataFrame
cpop.PLOT_VALUE = cpop.ZIP .|> x -> get(jdf.IDL_ZIPCODE |> countmap, x, 0)
# Downloaded from https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
counties = gmtread("/data/ursa_research/ris3/Ashlin/cb_2018_us_zcta510_500k/cb_2018_us_zcta510_500k.shp")
dfc = DataFrame(ZIP = map(x->x.attrib["ZCTA5CE10"],counties),ORDER=1:length(counties))
joineddata = @chain leftjoin(dfc,cpop,on= [:ZIP],makeunique=true) begin
	@orderby(:ORDER)
end
cptvallog = makecpt(range=(0,130),C=:matter)
# I think these are the options: 
# https://docs.juliaplots.org/latest/generated/colorschemes/#cmocean
# algae looks nice, but it doesn't show contrast well
# speed is okay
# I like deep
joineddata.PLOT_VALUE = replace(joineddata.PLOT_VALUE,missing => 0)
GMT.plot(counties,level=joineddata.PLOT_VALUE,cmap=cptvallog,close=true,fill="+z",pen=0.25,region=(-71.9,-71.1,41.1,42.05),
           proj=:guess,colorbar=true,figname="choroplethlog.png",title="Number of sequences by ZIP code")

#-----------------------------------------------------------
@info "Writing output file"
#-----------------------------------------------------------

CSV.write("out.csv", jdf, quotestrings=true)

#===========================================================
# Example commands
===========================================================#

RUN_EXAMPLES = false
if(RUN_EXAMPLES)
	# What are the column names in this DataFrame?
	names(vdf)

	# What are ALL the column names in this DataFrame?
	names(vdf) |> show

	# What are the possible values in this column?
	vdf.IDL_RACE |> unique
	ddf[!,"Collection Method"] |> unique # Useful if column name is not a valid symbol (e.g., contains a space)

	# What is the data dictionary information for the field :IDL_AGE?
		ddf[ddf[!,"Variable Name"] .== "IDL_AGE",:]

	# See most frequent values in jdf.IDL_ZIPCODE, sorted by number of occurrences
	first(sort(combine(groupby(jdf, :IDL_ZIPCODE), nrow), :nrow, rev=true), 15)

	# Select only females
	#filter(:IDL_SEX => ==("Female"), jdf) # This is standard, but it fails with missing values
	filter(:IDL_SEX => x -> !ismissing(x) && x == "Female", jdf)

	# Select by age
	filter(:IDL_AGE => x -> !ismissing(x) && 20 ≤ x < 30, jdf)

	# Simple plots
	vdf.IDL_AGE |> histogram
	jdf.IDL_AGE |> boxplot

	# Look at 10 random accession numbers
	sample(vdf.ACCESSION_NUMBER, 10)

	# See accession number with typo (extra dash) ...
	temp = filter(x -> !ismissing(x) && first(x,4) == "RI--", vdf.ACCESSION_NUMBER)
	# ... And see that it doesn't appear in the gisaid data
	filter(:ACCESSION_NUMBER => x -> !ismissing(x) && x == temp, gdf)

	# Compare sets of column names
	x = names(vdf);
	y = ddf[!, "Variable Name"];
	x ∩ y # Columns found in the data dictionary
	x[x.∉Ref(y)] # Columns not found in the data dictionary
	y[y.∉Ref(x)] # Dictionary entries not found in the variants columns

	# Verify that accession numbers map to a single unique row in the variant data
	unique(1:length(jdf.ACCESSION_NUMBER) .|> y->filter(:ACCESSION_NUMBER => x-> !ismissing(x) && x==jdf.ACCESSION_NUMBER[y], vdf) |> nrow) == [1]

	# Find any person IDs that occur more than once
	filter(x->last(x)>1, countmap(jdf.IDL_PERSON_ID))

end

