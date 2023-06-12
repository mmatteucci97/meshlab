# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


set(SOURCES filter_cut_and_fill.cpp)

set(HEADERS filter_cut_and_fill.h)

add_meshlab_plugin(filter_cut_and_fill ${SOURCES} ${HEADERS})

target_link_libraries(filter_cut_and_fill PRIVATE OpenGL::GLU)
