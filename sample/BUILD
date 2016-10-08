package(
    default_visibility = ["//visibility:public"],
)

cc_binary(
    name = "sample",
    srcs = ["sample.cc"],
    deps = [
        "//indexing/rct/public:dense_vec_data",
    ],
    linkopts = [
        "-lglog",
    ],
)
