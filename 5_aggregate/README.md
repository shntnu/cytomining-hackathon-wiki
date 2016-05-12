# Aggregating cell-level profiles to create treatment profiles

@mrohban
@Swarchal

We created a function that aggregates data to a user specified level, i.e from single cell data to well averrages. The method of aggregation is up to the user, though we have the default set as median. The `aggregate_data()` function works on in-memory dataframes and out-of-memory sql tables.
