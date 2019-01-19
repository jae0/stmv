stmv_error_codes = function() {
  stmv_error = list()
  stmv_error["todo"] = 0L
  stmv_error["complete"] = 1L
  stmv_error["outside_bounds"] = 2L
  stmv_error["too_shallow"] = 3L
  stmv_error["prediction_area"] = 4L
  stmv_error["insufficient_data"] = 5L
  stmv_error["variogram_failure"] = 6L
  stmv_error["variogram_range_limit"] = 7L
  stmv_error["prediction_error"] = 8L
  stmv_error["prediction_update_error"] = 9L
  stmv_error["statistics_update_error"] = 10L
  stmv_error["local_model_error"] = 11L
  stmv_error["unknown"] = -1L
  return(stmv_error)
}
