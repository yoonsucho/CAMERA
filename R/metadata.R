CAMERA$set("public", "get_metadata", function(exposure_ids=self$exposure_ids, outcome_ids=self$outcome_ids) {
    self$exposure_metadata <- ieugwasr::gwasinfo(exposure_ids)
    self$outcome_metadata <- ieugwasr::gwasinfo(outcome_ids)
    return(list(self$exposure_metadata, self$outcome_metadata))
})
