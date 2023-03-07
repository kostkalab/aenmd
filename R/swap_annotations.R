
#- here we are swapping the ad_* function
#  aenmd uses to access annotation data.

#- since they are exported, we need to swap 
#  then in two environments namespace:aenmd 
#  and package:aenmd.
#  see: https://cran.r-project.org/doc/manuals/R-ints.html#Namespaces

#- Note: this works, but there should be a better way.

swap_fu <- function(fname, PACKAGENAME){
#---------------------------------------

    #- this takes care of the function in aenmd namespace
    #  i.e., environment namespace:aenmd
    rlang::env_binding_unlock(env = asNamespace('aenmd'),
                              nms = fname)
    assign(x     = fname, 
           value = get0(fname, envir = asNamespace(PACKAGENAME)),
           envir = asNamespace('aenmd'))
 
    rlang::env_binding_lock(env = asNamespace('aenmd'),
                            nms = fname)

    #- this takes care of the package:aenmd environment
    #  where the exported functions have a second copy
    #  (why are they copied and not pointed to??)
    aenmdenv <- rlang::search_env('package:aenmd')
    rlang::env_binding_unlock(env = aenmdenv,
                              nms = fname)
    assign(x     = fname, 
           value = get0(fname, envir = asNamespace(PACKAGENAME)),
           envir = aenmdenv)

    rlang::env_binding_lock(env = aenmdenv,
                            nms = fname)
}

#' Swap out annotation package used by aenmd
#' @param PACKAGENAME String. Name of the new annotation package to be used.
#' @return None. This function is called only for the side effect of swapping
#' out the annotation package aenmd is using.
#' @export
ad_swap_annotation <- function(PACKAGENAME){
#--------------------------------------------

    #- load namespace of new annotation package
    requireNamespace(PACKAGENAME, quietly = TRUE)

    #- swap out aenmd's ad_* functions
    swap_fu('ad_get_genome',       PACKAGENAME)
    swap_fu('ad_get_txs',          PACKAGENAME)
    swap_fu('ad_get_txs_mask',     PACKAGENAME)
    swap_fu('ad_get_spl_mask',     PACKAGENAME)
    swap_fu('ad_get_cds_by_tx',    PACKAGENAME)
    swap_fu('ad_is_single_exn_tx', PACKAGENAME)
    swap_fu('ad_is_ptc_snv',       PACKAGENAME)

    #- adjust the data package name
    old_package <- ._EA_dataPackage_name
    rlang::env_binding_unlock(env = asNamespace('aenmd'),
                              nms = '._EA_dataPackage_name')
    assign(x     = '._EA_dataPackage_name', 
           value = PACKAGENAME,
           envir = asNamespace('aenmd'))
 
    rlang::env_binding_lock(env = asNamespace('aenmd'),
                            nms = '._EA_dataPackage_name')

    #- unload the old package
    unloadNamespace(asNamespace(old_package))
}


