#' Return column names matching a prefix
#'
#' Useful for returning names of metadata or feature data columns in a given
#' dataframe or database table
#'
#' @param x dataframe o database table
#' @param prefix string, prefix of columns
#'
#' @examples
#' data(iris)
#' # get the column names  in iris beginning with "S"
#' get_cols(iris, prefix = "S")
#' @export

get_cols <- function(x, prefix) {
    colnames(x)[grep(sprintf("^%s", prefix), colnames(x))]
}


#' aggregate data from single cell measurements
#' 
#' This functionn aggregates single cell level data to a user-specified level,
#' typically to an image or well median for each feature.
#' 
#' @param data dataframe or sql table
#' @param grouping.column string, name of column with which to group by
#' @param feature.columns strings, names of featuredata columns
#' @param meta.columns strings, names of metadata columns
#' @param aggregation.function function with which to aggregate data
#' 
#' @import dplyr
#' 
#' @examples
#' # add example
#' # ideally need some example dataset
#' 
#' @export

aggregate_data <- function(data, grouping.column, feature.columns, meta.columns,
			   aggregation.function = median) {
    
    # compute aggregates of the feature data
    data.aggr <- data %>%
	group_by_(grouping.column) %>%
	select(one_of(feature.columns)) %>%
	summarise_each(funs(aggregation.function))
    
    # collapses specified metadata to the same level as feature data
    meta.data <- data %>%
	select(one_of(meta.columns)) %>%
	unique()
    
    # merge featuredata and metadata by the grouping column
    # if performed on a database will collect queries, if performed on a
    # dataframe will return the dataframe
    inner_join(x = data.aggr, y = meta.data, by = grouping.column) %>%
	collect()
}
