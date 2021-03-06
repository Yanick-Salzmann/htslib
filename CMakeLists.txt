cmake_minimum_required(VERSION 3.13)
project(htslib)

include(ExternalProject)

set(ZLIB_ROOT ${CMAKE_BINARY_DIR}/zlib)

ExternalProject_Add(
        zlib
        PREFIX ${ZLIB_ROOT}
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG master
        UPDATE_COMMAND ""
        PATCH_COMMAND ""
        SOURCE_DIR ${CMAKE_BINARY_DIR}/src/zlib
        INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/zlib
        CMAKE_CACHE_ARGS
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/zlib
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
)

add_library(htslib SHARED
        bzip2/blocksort.c
        bzip2/bzlib.c
        bzip2/bzlib.h
        bzip2/bzlib_private.h
        bzip2/compress.c
        bzip2/crctable.c
        bzip2/decompress.c
        bzip2/huffman.c
        bzip2/randtable.c
        cram/cram.h
        cram/cram_codecs.c
        cram/cram_codecs.h
        cram/cram_decode.c
        cram/cram_encode.c
        cram/cram_decode.h
        cram/cram_encode.h
        cram/cram_external.c
        cram/cram_index.c
        cram/cram_index.h
        cram/cram_io.c
        cram/cram_io.h
        cram/cram_samtools.c
        cram/cram_samtools.h
        cram/cram_stats.c
        cram/cram_stats.h
        cram/mFILE.c
        cram/mFILE.h
        cram/misc.h
        cram/open_trace_file.c
        cram/open_trace_file.h
        cram/os.h
        cram/pooled_alloc.c
        cram/pooled_alloc.h
        cram/rANS_byte.h
        cram/rANS_static.c
        cram/rANS_static.h
        cram/string_alloc.c
        cram/string_alloc.h
        config.h
        bcf_sr_sort.c
        bcf_sr_sort.h
        bgzf.c
        bgzip.c
        errmod.c
        faidx.c
        header.c
        header.h
        hfile.c
        hfile_gcs.c
        hfile_internal.h
        #hfile_libcurl.c
        #hfile_net.c
        #hfile_s3.c
        #hfile_s3_write.c
        hts.c
        hts_internal.h
        hts_os.c
        htsfile.c
        kfunc.c
        #knetfile.c
        kstring.c
        md5.c
        multipart.c
        #plugin.c
        probaln.c
        realn.c
        regidx.c
        region.c
        sam.c
        sam_internal.h
        synced_bcf_reader.c
        tabix.c
        tbx.c
        textutils.c
        textutils_internal.h
        thread_pool.c
        thread_pool_internal.h
        vcf.c
        vcf_sweep.c
        vcfutils.c

        liblzma/api/lzma/base.h
        liblzma/api/lzma/bcj.h
        liblzma/api/lzma/block.h
        liblzma/api/lzma/check.h
        liblzma/api/lzma/container.h
        liblzma/api/lzma/delta.h
        liblzma/api/lzma/filter.h
        liblzma/api/lzma/hardware.h
        liblzma/api/lzma/index.h
        liblzma/api/lzma/index_hash.h
        liblzma/api/lzma/lzma.h
        liblzma/api/lzma/stream_flags.h
        liblzma/api/lzma/version.h
        liblzma/api/lzma/vli.h
        liblzma/api/lzma.h
        liblzma/check/check.h
        liblzma/check/check.c
        liblzma/check/crc32_fast.c
        liblzma/check/crc32_table.c
        liblzma/check/crc32_table_be.h
        liblzma/check/crc32_table_le.h
        liblzma/check/crc64_fast.c
        liblzma/check/crc64_table.c
        liblzma/check/crc64_table_be.h
        liblzma/check/crc64_table_le.h
        liblzma/check/crc_macros.h
        liblzma/check/sha256.c
        liblzma/common/alone_decoder.c
        liblzma/common/alone_decoder.h
        liblzma/common/alone_encoder.c
        liblzma/common/auto_decoder.c
        liblzma/common/block_buffer_decoder.c
        liblzma/common/block_buffer_encoder.c
        liblzma/common/block_decoder.c
        liblzma/common/block_decoder.h
        liblzma/common/block_encoder.c
        liblzma/common/block_encoder.h
        liblzma/common/block_header_decoder.c
        liblzma/common/block_header_encoder.c
        liblzma/common/block_util.c
        liblzma/common/common.c
        liblzma/common/common.h
        liblzma/common/easy_buffer_encoder.c
        liblzma/common/easy_decoder_memusage.c
        liblzma/common/easy_encoder.c
        liblzma/common/easy_encoder_memusage.c
        liblzma/common/easy_preset.c
        liblzma/common/easy_preset.h
        liblzma/common/filter_buffer_decoder.c
        liblzma/common/filter_buffer_encoder.c
        liblzma/common/filter_common.c
        liblzma/common/filter_common.h
        liblzma/common/filter_decoder.c
        liblzma/common/filter_decoder.h
        liblzma/common/filter_encoder.c
        liblzma/common/filter_encoder.h
        liblzma/common/filter_flags_decoder.c
        liblzma/common/filter_flags_encoder.c
        liblzma/common/hardware_physmem.c
        liblzma/common/index.h
        liblzma/common/index.c
        liblzma/common/index_decoder.c
        liblzma/common/index_encoder.c
        liblzma/common/index_hash.c
        liblzma/common/mythread.h
        liblzma/common/stream_buffer_decoder.c
        liblzma/common/stream_buffer_encoder.c
        liblzma/common/stream_decoder.c
        liblzma/common/stream_decoder.h
        liblzma/common/stream_encoder.c
        liblzma/common/stream_encoder.h
        liblzma/common/stream_flags_common.c
        liblzma/common/stream_flags_common.h
        liblzma/common/stream_flags_decoder.c
        liblzma/common/stream_flags_encoder.c
        liblzma/common/sysdefs.h
        liblzma/common/tuklib_common.h
        liblzma/common/tuklib_config.h
        liblzma/common/tuklib_cpucores.c
        liblzma/common/tuklib_cpucores.h
        liblzma/common/tuklib_exit.c
        liblzma/common/tuklib_exit.h
        liblzma/common/tuklib_gettext.h
        liblzma/common/tuklib_integer.h
        liblzma/common/tuklib_mbstr.h
        liblzma/common/tuklib_mbstr_fw.c
        liblzma/common/tuklib_mbstr_width.c
        liblzma/common/tuklib_open_stdxxx.c
        liblzma/common/tuklib_open_stdxxx.h
        liblzma/common/tuklib_physmem.c
        liblzma/common/tuklib_physmem.h
        liblzma/common/tuklib_progname.c
        liblzma/common/tuklib_progname.h
        liblzma/common/vli_decoder.c
        liblzma/common/vli_encoder.c
        liblzma/common/vli_size.c
        liblzma/delta/delta_common.c
        liblzma/delta/delta_common.h
        liblzma/delta/delta_decoder.c
        liblzma/delta/delta_decoder.h
        liblzma/delta/delta_encoder.c
        liblzma/delta/delta_encoder.h
        liblzma/delta/delta_private.h
        liblzma/lz/lz_decoder.c
        liblzma/lz/lz_decoder.h
        liblzma/lz/lz_encoder.c
        liblzma/lz/lz_encoder.h
        liblzma/lz/lz_encoder_hash.h
        liblzma/lz/lz_encoder_hash_table.h
        liblzma/lz/lz_encoder_mf.c
        liblzma/lzma/fastpos.h
        liblzma/lzma/fastpos_table.c
        liblzma/lzma/lzma2_decoder.c
        liblzma/lzma/lzma2_decoder.h
        liblzma/lzma/lzma2_encoder.c
        liblzma/lzma/lzma2_encoder.h
        liblzma/lzma/lzma_common.h
        liblzma/lzma/lzma_decoder.c
        liblzma/lzma/lzma_decoder.h
        liblzma/lzma/lzma_encoder.c
        liblzma/lzma/lzma_encoder.h
        liblzma/lzma/lzma_encoder_optimum_fast.c
        liblzma/lzma/lzma_encoder_optimum_normal.c
        liblzma/lzma/lzma_encoder_presets.c
        liblzma/lzma/lzma_encoder_private.h
        liblzma/rangecoder/price.h
        liblzma/rangecoder/price_table.c
        liblzma/rangecoder/price_table.c
        liblzma/rangecoder/range_common.h
        liblzma/rangecoder/range_decoder.h
        liblzma/rangecoder/range_encoder.h
        liblzma/simple/arm.c
        liblzma/simple/armthumb.c
        liblzma/simple/ia64.c
        liblzma/simple/powerpc.c
        liblzma/simple/simple_coder.c
        liblzma/simple/simple_coder.h
        liblzma/simple/simple_decoder.c
        liblzma/simple/simple_decoder.h
        liblzma/simple/simple_encoder.c
        liblzma/simple/simple_encoder.h
        liblzma/simple/simple_private.h
        liblzma/simple/sparc.c
        liblzma/simple/x86.c

        pthread/src/barrier.c
        pthread/src/barrier.h
        pthread/src/clock.c
        pthread/src/cond.c
        pthread/src/cond.h
        pthread/src/misc.c
        pthread/src/misc.h
        pthread/src/mutex.c
        pthread/src/mutex.h
        pthread/src/nanosleep.c
        pthread/src/ref.c
        pthread/src/ref.h
        pthread/src/rwlock.c
        pthread/src/rwlock.h
        pthread/src/sched.c
        pthread/src/sem.c
        pthread/src/sem.h
        pthread/src/spinlock.c
        pthread/src/thread.c
        pthread/src/thread.h
        pthread/src/winpthread_internal.h
        pthread/src/wpth_ver.h
        )

add_dependencies(htslib zlib)

target_compile_definitions(htslib PUBLIC -DHTS_BUILDING_LIBRARY)
target_include_directories(htslib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/zlib/include)
target_link_libraries(htslib ${ZLIB_LIBRARY} ${CMAKE_CURRENT_BINARY_DIR}/zlib/lib/libzlibstatic.a ws2_32.lib)