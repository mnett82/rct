package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "public",
    deps = [
        ":dist_data",
        ":rct",
    ],
)

cc_library(
    name = "dist_data",
    hdrs = ["dist_data.h"],
)

cc_library(
    name = "rct",
    srcs = ["rct.cc"],
    hdrs = ["rct.h"],
    deps = [
        ":dist_data",
    ],
    linkopts = [
        "-lglog",
    ],
)
