
# ------------------------------------------------------------------
# prepare raw deaths linelist data
#' @export

prepare_deaths <- function(deaths) {
    
    # subset to NHSE deaths only
    deaths <- subset(deaths, nhse == "reported by NHSE")
    
    # make age numeric
    deaths$age <- as.numeric(as.character(deaths$age))
    
    # make sex boolean
    deaths$sex_female <- c(TRUE, FALSE)[match(deaths$gender, c("female", "male"))]
    
    # add dummy region column if not present
    if (!("region" %in% names(deaths))) {
        deaths$region <- NA
    }
    
    # subset and rename columns
    deaths <- data.frame(age = deaths$age,
                         sex_female = deaths$sex_female,
                         region = deaths$region,
                         date_onset = deaths$onsetdate,
                         date_admission = deaths$dateadmission,
                         date_swab = deaths$dateswab_nhse,
                         date_labtest = deaths$dateofresult_nhse,
                         date_death = deaths$dod)
    
    return(deaths)
}

# ------------------------------------------------------------------
# prepare raw SitRep data
#' @export

prepare_sitrep <- function(sitrep, deaths) {
    
    # sum sitrep over all locations
    sitrep$all <- rowSums(subset(sitrep, select = -c(date, metric_name)))
    
    # get into wide format
    sitrep_wide <- tidyr::pivot_wider(sitrep, id_cols = date,
                                      names_from = metric_name, values_from = all)
    
    # subset and rename columns
    sitrep <- data.frame(
        date = sitrep_wide$date,
        total_hdu = sitrep_wide$number_of_confirmed_covid_19_patients_in_hdu_at_0800_total,
        total_itu = sitrep_wide$number_of_confirmed_covid_19_patients_in_itu_at_0800_total,
        total_hdu_icu = sitrep_wide$number_of_confirmed_covid_19_patients_in_hdu_itu_at_0800_total,
        total_idu = sitrep_wide$number_of_confirmed_covid_19_patients_in_infectious_disease_unit_beds_at_0800_total,
        total_other = sitrep_wide$number_of_confirmed_covid_19_patients_in_any_other_beds_at_0800_total,
        new_discharges = sitrep_wide$number_of_covid_19_discharges_in_the_past_24_hours_total,
        new_discharges_usual = sitrep_wide$number_of_covid_19_discharges_to_usual_place_of_residence_in_the_last_24_hours_total,
        new_inpatients_diagnosed = sitrep_wide$number_of_inpatients_diagnosed_with_covid_19_in_last_24_hours_total,
        new_admissions = sitrep_wide$number_of_patients_admitted_with_covid_19_in_last_24_hours_total,
        swab_waiting = sitrep_wide$number_of_patients_that_have_had_diagnostic_swabbing_for_covid_19_and_are_awaiting_results_at_0800)
    
    # create column for total general beds occupied
    sitrep$total_general = sitrep$total_idu + sitrep$total_other
    
    #### combine with deaths data
    # move death dates forward by one day to correspond with SitRep
    # (which is recorded at 8:00am for previous 24h period)
    deaths$date_death <- deaths$date_death + 1
    
    # add deaths column to processed SitRep
    tab1 <- table(deaths$date_death)
    sitrep$deaths <- as.vector(tab1)[match(sitrep$date, as.Date(names(tab1)))]
    
    # create column for numeric date
    sitrep$date_numeric <- as.numeric(sitrep$date - min(sitrep$date))
    
    return(sitrep)
}

# ------------------------------------------------------------------
# prepare CHESS data
#' @export

prepare_chess <- function(chess) {
    
    # make sex boolean
    chess$sex_female <- c(TRUE, FALSE)[match(chess$sex, c("female", "male"))]
    
    # re-level final outcome. Consider transfers equivalent to discharges
    chess$final_outcome <- c(c("death", "discharge", "discharge"))[match(chess$finaloutcome, c("Death", "Discharged", "Transfered"))]
    
    # create dummy column for date updated
    chess$date_updated = as.Date(NA)
    
    # subset and rename columns
    chess <- data.frame(case_id = chess$caseid,
                        age = chess$ageyear,
                        sex_female = chess$sex_female,
                        date_onset = chess$estimateddate,
                        date_admission = chess$hospitaladmissiondate,
                        date_swab = chess$infectionswabdate,
                        date_labtest = chess$labtestdate,
                        date_icu = chess$dateadmittedicu,
                        date_leave_icu = chess$dateleavingicu,
                        date_final_outcome = chess$finaloutcomedate,
                        final_outcome = chess$final_outcome,
                        date_updated = chess$date_updated)
    
    return(chess)
}

# ------------------------------------------------------------------
# prepare COCIN data
#' @export

prepare_cocin <- function(cocin) {
    
    # create dummy labtest column
    if (!("date_labtest" %in% names(cocin))) {
        cocin$date_labtest <- as.Date(NA)
    }
    
    # create dummy swab column
    if (!("date_swab" %in% names(cocin))) {
        cocin$date_swab <- as.Date(NA)
    }
    
    # fix date_last_assess format
    #cocin$date_last_assess <- as.Date(cocin$date_last_assess, origin = as.Date("1970-01-01"))
    
    # fix other date formats
    cols <- c("date_onset", "date_admission", "date_swab", "date_icu",
              "date_leave_icu", "date_final_outcome", "date_last_assess")
    for (i in seq_along(cols)) {
        cocin[,cols[i]] <- as.Date(as.character( cocin[,cols[i]] ))
    }
    
    # subset and rename columns
    cocin <- data.frame(case_id = cocin$case_id,
                        age = cocin$age_year,
                        sex_female = cocin$sex_female,
                        date_onset = cocin$date_onset,
                        date_admission = cocin$date_admission,
                        date_swab = cocin$date_swab,
                        date_labtest = cocin$date_labtest,
                        date_icu = cocin$date_icu,
                        date_leave_icu = cocin$date_leave_icu,
                        date_final_outcome = cocin$date_final_outcome,
                        final_outcome = cocin$final_outcome,
                        date_updated = cocin$date_last_assess)
    
    return(cocin)
}

# ------------------------------------------------------------------
# prepare COCIN data sent by Lucy (with regions)
#' @export

prepare_cocin_v2 <- function(cocin) {
    
    # create dummy labtest column
    cocin$date_swab <- as.Date(NA)
    cocin$date_labtest <- as.Date(NA)
    
    # fix other date formats
    cols <- c("date_onset", "date_admission", "date_swab", "date_icu",
              "date_leave_icu", "date_final_outcome", "date_last_assess")
    for (i in seq_along(cols)) {
        cocin[,cols[i]] <- as.Date(as.character( cocin[,cols[i]] ))
    }
    
    # subset and rename columns
    cocin <- data.frame(case_id = cocin$case_id,
                        age = cocin$age_year,
                        sex_female = cocin$sex_female,
                        region = tolower(cocin$region),
                        date_onset = cocin$date_onset,
                        date_admission = cocin$date_admission,
                        date_swab = as.Date(NA),#cocin$date_swab,
                        date_labtest = cocin$date_labtest,
                        date_icu = cocin$date_icu,
                        date_leave_icu = cocin$date_leave_icu,
                        date_final_outcome = cocin$date_final_outcome,
                        final_outcome = cocin$final_outcome,
                        date_updated = cocin$date_last_assess)
    
    return(cocin)
}

# ------------------------------------------------------------------
# prepare individual-level data
#' @export

prepare_indlevel <- function(indlevel) {
    
    # remove duplicate IDs
    # NB. these cases can represent multiple trips to ICU
    dup_ID <- indlevel$case_id[duplicated(indlevel$case_id)]
    indlevel <- subset(indlevel, !(indlevel$case_id %in% dup_ID))
    
    # replace labtest date with admission date if labest < admission
    w <- which(indlevel$date_labtest < indlevel$date_admission)
    indlevel$date_labtest[w] <- indlevel$date_admission[w]
    
    # if have final outcome, must have accompanying date
    indlevel <- subset(indlevel, !(!is.na(final_outcome) & is.na(date_final_outcome)))
    
    # record whether went via ICU
    indlevel$icu <- !is.na(indlevel$date_icu)
    
    # get date stepped down from ICU. Assume that discharges from ICU without a
    # specified date of leaving ICU have stepdown for 0 days.
    indlevel$date_stepdown <- as.Date(NA)
    # discharged from ICU without date leaving ICU. Assume stepdown on date of final outcome
    w <- which(indlevel$icu &
               is.na(indlevel$date_leave_icu) & 
               (indlevel$final_outcome == "discharge"))
    indlevel$date_stepdown[w] <- indlevel$date_final_outcome[w]
    # discharged from ICU with date leaving ICU. Assume stepdown on date of leaving ICU
    w <- which(indlevel$icu &
               !is.na(indlevel$date_leave_icu) &
               (indlevel$final_outcome == "discharge"))
    indlevel$date_stepdown[w] <- indlevel$date_leave_icu[w]
    # registered date of leaving ICU, and date of death not equal to date leaving ICU. Assume stepdown on date of leaving ICU
    w <- which(indlevel$icu &
               !is.na(indlevel$date_leave_icu) &
               (indlevel$final_outcome == "death") &
               (indlevel$date_leave_icu != indlevel$date_final_outcome))
    indlevel$date_stepdown[w] <- indlevel$date_leave_icu[w]
    
    # record whether stepdown occurred
    indlevel$stepdown <- NA
    indlevel$stepdown[indlevel$icu] <- !is.na(indlevel$date_stepdown[indlevel$icu])
    
    # filter out deaths from stepdown
    #table(indlevel$final_outcome[which(indlevel$stepdown == TRUE)])
    indlevel <- subset(indlevel, !( !is.na(stepdown) & (stepdown == TRUE) & (final_outcome == "death") ))
    
    # ensure dates in correct sequence
    indlevel <- subset(indlevel, !icu | date_admission <= date_icu)
    indlevel <- subset(indlevel, is.na(stepdown) | !stepdown | date_icu <= date_stepdown)
    indlevel <- subset(indlevel, is.na(date_final_outcome) | date_admission <= date_final_outcome)
    indlevel <- subset(indlevel, !icu | is.na(date_final_outcome) | date_icu <= date_final_outcome)
    indlevel <- subset(indlevel, is.na(stepdown) | !stepdown | is.na(date_final_outcome) | date_stepdown <= date_final_outcome)
    
    # add censoring date
    indlevel$date_censor <- as.Date(NA)
    for (i in seq_len(nrow(indlevel))) {
        #tmp <- c(indlevel$date_admission[i], indlevel$date_swab[i], indlevel$date_labtest[i],
        #         indlevel$date_icu[i], indlevel$date_leave_icu[i], indlevel$date_final_outcome[i],
        #         indlevel$date_updated[i])
        tmp <- indlevel$date_updated[i]
        if (any(!is.na(tmp))) {
            indlevel$date_censor[i] <- max(tmp, na.rm = TRUE)
        }
    }
    
    return(indlevel)
}

# ------------------------------------------------------------------
# prepare raw SitRep data in age-stratified format
#' @export

prepare_sitrep_age <- function(sitrep, deaths) {
    
    # define age bands
    age_bands <- data.frame(group = 1:5,
                            lower = c(0, 6, 18, 65, 85),
                            name_sitrep = c("_0_5", "_6_17", "_18_64", "_65_84", "_85"),
                            stringsAsFactors = FALSE)
    
    # sum sitrep over all locations
    sitrep$all <- rowSums(subset(sitrep, select = -c(date, metric_name)))
    
    # subset metric names
    keep_names <- c("number_of_confirmed_covid_19_patients_in_hdu_itu_at_0800",
                    "number_of_confirmed_covid_19_patients_in_infectious_disease_unit_beds_at_0800",
                    "number_of_confirmed_covid_19_patients_in_any_other_beds_at_0800",
                    "number_of_covid_19_discharges_in_the_past_24_hours",
                    "number_of_inpatients_diagnosed_with_covid_19_in_last_24_hours",
                    "number_of_patients_admitted_with_covid_19_in_last_24_hours")
    
    sitrep <- sitrep[grep(paste(keep_names, collapse = "|"), sitrep$metric_name),]
    
    # get age group from metric name and drop suffix
    sitrep$age_group <- NA
    for (i in seq_len(nrow(age_bands))) {
        w <- grep(age_bands$name_sitrep[i], as.character(sitrep$metric_name))
        sitrep$metric_name[w] <- stringr::str_remove(sitrep$metric_name[w], age_bands$name_sitrep[i])
        sitrep$age_group[w] <- age_bands$group[i]
    }
    
    # only keep metric names with an associated age group (drop total counts)
    sitrep <- subset(sitrep, !is.na(age_group))
    
    # get into wide format
    sitrep_wide <- tidyr::pivot_wider(sitrep, id_cols = c(date, age_group), names_from = metric_name, values_from = all)
    
    # re-order rows and columns
    sitrep_wide <- sitrep_wide[order(sitrep_wide$date, sitrep_wide$age_group),]
    sitrep_wide <- subset(sitrep_wide, select = c("date", "age_group", keep_names))
    
    # rename columns
    new_names <- c("date",
                   "age_group",
                   "total_hdu_icu",
                   "total_idu",
                   "total_other",
                   "new_discharges",
                   "new_inpatients_diagnosed",
                   "new_admissions")
    names(sitrep_wide) <- new_names
    
    # create column for total general beds occupied
    sitrep_wide$total_general <- sitrep_wide$total_idu + sitrep_wide$total_other
    
    
    #### combine with deaths data
    # move death dates forward by one day to correspond with SitRep
    # (which is recorded at 8:00am for previous 24h period)
    deaths$date_death <- deaths$date_death + 1
    
    # get death age groups
    deaths$age_group <- as.numeric(cut(deaths$age, breaks = c(age_bands$lower, Inf), right = FALSE))
    
    # tabulate deaths by age and date
    tab1 <- as.data.frame(table(deaths$date_death, deaths$age_group))
    names(tab1) <- c("date", "age_group", "deaths")
    
    # merge with SitRep
    sitrep_wide <- merge(sitrep_wide, tab1, all.x = TRUE)
    
    # re-order rows and columns
    new_names <- c(new_names, "total_general", "deaths")
    sitrep_wide <- sitrep_wide[order(sitrep_wide$date, sitrep_wide$age_group),]
    sitrep_wide <- subset(sitrep_wide, select = new_names)
    
    # create column for numeric date
    sitrep_wide$date_numeric <- as.numeric(sitrep_wide$date - min(sitrep_wide$date))
    
    # create gaps in data for missing dates
    date_buffer <- expand.grid(date_numeric = 0:max(sitrep_wide$date_numeric),
                               age_group = seq_len(nrow(age_bands)))
    sitrep_wide <- merge(sitrep_wide, date_buffer, all = TRUE)
    
    # re-order rows and columns
    new_names <- c("date_numeric", new_names)
    sitrep_wide <- sitrep_wide[order(sitrep_wide$date_numeric, sitrep_wide$age_group),]
    sitrep_wide <- subset(sitrep_wide, select = new_names)
    
    # ensure no missing dates
    sitrep_wide$date <- as.Date(sitrep_wide$date_numeric, origin = as.Date(min(sitrep_wide$date, na.rm = TRUE)))
    
    # replace death NAs with 0
    sitrep_wide$deaths[is.na(sitrep_wide$deaths)] <- 0
    
    return(sitrep_wide)
}

