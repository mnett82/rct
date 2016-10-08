package(
    default_visibility = ["//visibility:public"],
)

cc_test(
    name = "dense_feature_vector_test",
    srcs = ["dense_feature_vector_test.cc"],
    size = "small",
    deps = [
        "//indexing/rct/public:dense_feature_vector",
        "//third_party/gtest:gtest_main",
    ],
)

cc_test(
    name = "metrics_test",
    srcs = ["metrics_test.cc"],
    size = "small",
    deps = [
        "//indexing/rct/public:dense_feature_vector",
        "//indexing/rct/public:metrics",
        "//third_party/gtest:gtest_main",
    ],
)
