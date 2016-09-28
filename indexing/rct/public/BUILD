package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "public",
    deps = [
        ":dist_data",
        ":random",  # TODO: this dependency needs to go.
        ":rct",
    ],
)

cc_library(
    name = "dist_data",
    hdrs = ["dist_data.h"],
)

# TODO: use standard PRNG instead.
cc_library(
    name = "random",
    srcs = ["random.cc"],
    hdrs = ["random.h"],
)

cc_library(
    name = "rct",
    srcs = ["rct.cc"],
    hdrs = ["rct.h"],
    deps = [
        ":dist_data",
        ":random",  # TODO: this dependency needs to go.
    ],
)
