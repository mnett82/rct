package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "dense_feature_vector",
    hdrs = ["dense_feature_vector.h"],
)

cc_library(
    name = "metrics",
    hdrs = ["metrics.h"],
    deps = [
        ":dense_feature_vector",
    ],
)
