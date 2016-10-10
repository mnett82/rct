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

cc_library(
    name = "rank_cover_tree",
    srcs = ["rank_cover_tree-inl.h"],
    hdrs = ["rank_cover_tree.h"],
    # TODO: Bring in GLog as a source dependency.
    linkopts = [
        "-lglog",
    ],
)
