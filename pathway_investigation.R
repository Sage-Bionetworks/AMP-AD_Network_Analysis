
geneRef <- unique(unlist(dbgap))

foo_foo<- utilityFunctions::outerSapplyParallel(utilityFunctions::fisherWrapperPval,
                                                dbgap,
                                                dbgap,
                                                geneRef)

foo_foo2<- utilityFunctions::outerSapplyParallel(utilityFunctions::fisherWrapperOR,
                                                dbgap,
                                                dbgap,
                                                geneRef)
