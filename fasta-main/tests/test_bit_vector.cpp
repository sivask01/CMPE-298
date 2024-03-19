/**
 * Copyright (c) David C. Anastasiu & Gheorghi Guzun
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */
#include <iostream>     // std::cout
#include <limits.h>
#include <chrono>
#include <vector>
#include <gtest/gtest.h>

//header file
#include <fasta/utils/BitVector.h>


using namespace fasta;

/**
* Test if vector is allocated with all 0s
*/
TEST(BitVector, zero_vector){
    uint64_t size = 10000;
    BitVector bv = BitVector(size);
    EXPECT_EQ(bv.allocation_size(size),((size + 63) / 64));
    EXPECT_EQ(bv.isset(0),0);
    EXPECT_EQ(bv.isset(size-1),0);
    EXPECT_EQ(bv[0], 0);
}
/**
* Test various setting functions of vector
*/
TEST(BitVector, set_vector){
    uint64_t size = 10000;
    BitVector bv = BitVector(size);
    //set to 1
    for(uint64_t i = 1; i < size; i <<= 2){
        bv.set(i);
        EXPECT_TRUE(bv.isset(i));
    }
    //set some bits to 0
    for(uint64_t i = 1; i < size; i <<=4){
        bv.zero(i);
        EXPECT_FALSE(bv.isset(i));
    }
    for(uint64_t i = 1; i < size; i <<= 2){
        bv.set(i);
        EXPECT_TRUE(bv.isset(i));
    }
    bv.reset(); // reset all bits to 0
    EXPECT_FALSE(bv.isset(0));
    EXPECT_FALSE(bv.isset(size-1));
    EXPECT_FALSE(bv.isset(1024));

    BitVector bv2 = BitVector(size, true);
    for(uint64_t i = 1; i < size; i <<= 2){
        EXPECT_TRUE(bv2.isset(i));
    }
    //set some bits to 0
    for(uint64_t i = 1; i < size; i <<=4){
        bv2.zero(i);
        EXPECT_FALSE(bv2.isset(i));
    }
    bv2.reset(); // reset all bits to 1
    for(uint64_t i = 1; i < size; i <<= 2){
        EXPECT_TRUE(bv2.isset(i));
    }
    //set some bits to 0 and then set the blocks to 1
    for(uint64_t i = 1; i < size; i <<=4){
        bv2.zero(i);
        EXPECT_FALSE(bv2.isset(i));
        bv2.set_all(i);
        for(uint64_t j=i-i%64; j < i-i%64+64; ++j){
            EXPECT_TRUE(bv2.isset(j));
        }
    }
    //set some bits to 1 and then set the blocks to 0
    for(uint64_t i = 1; i < size; i <<=4){
        bv2.set(i);
        EXPECT_TRUE(bv2.isset(i));
        bv2.zero_all(i);
        for(uint64_t j=i-i%64; j < i-i%64+64; ++j){
            EXPECT_FALSE(bv2.isset(j));
        }
    }

    BitVector bv3 = BitVector(size);
    std::vector<uint32_t> elems;
    for(uint32_t i = 1; i < size; i <<=4){
        bv3.set(i);
        EXPECT_TRUE(bv3.isset(i));
        elems.push_back(i);
    }
    bv3.reset_all(elems.size(), elems.data());
    for(uint64_t i = 1; i < size; i <<=4){
        if(bv3.isset(i)){
            std::cout << i << std::endl;
        }
        EXPECT_FALSE(bv3.isset(i));
    }
}

/**
* Test vector index operator overload.
* Check if index operator correctly reads and writes values.
*/
TEST(BitVector, index_operator){
    uint64_t size = 128; //sz = 2
    BitVector bv = BitVector(size);
    bv.set(63);
    bv.set(62);
    bv.set(0);
    //Expect index operator to return the value of first index
    EXPECT_EQ(bv[0],3+((uint64_t)1<<63));

    //change value of first index
    bv[0] = 0;
    EXPECT_EQ(bv[0],0);
    EXPECT_EQ(bv.isset(1),0);

}

/**
* Test vector functionality with out of range inout
*/
TEST(BitVector, out_of_range_input){
    BitVector bv = BitVector(128);
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wshift-count-overflow"
    bv[0] = (uint64_t)1<<65;
    #pragma GCC diagnostic pop
    //truncated
    EXPECT_EQ(bv[0],0);
}

/**
* Test vector functionality with negative value. Also tests edge with max value
* Check that values are properly changed from negative int to uint64
*/
TEST(BitVector, negative_input){
    BitVector bv = BitVector(128);
    //setting to -1 should cause bv[0] to be max int of 64 bit int
    bv[0] = -1;
    EXPECT_EQ(bv.isset(0),1);
    EXPECT_EQ(bv.isset(63),1);
    EXPECT_EQ(bv.isset(1),1);

    // max 64bit uint = ULONG_MAX
    EXPECT_EQ(bv[0], ULONG_MAX);

    bv[1] = -2;
    EXPECT_EQ(bv[1], ULONG_MAX-1);
    EXPECT_EQ(bv[0], ULONG_MAX-1);

    bv[64] = -10;
    EXPECT_EQ(bv[64], ULONG_MAX-9);
    EXPECT_EQ(bv[0], ULONG_MAX-1);
}

/**
 * Test constructors
 */

TEST(BitVector, true_constructor){
    BitVector * bv = new BitVector(128, true);
    for(uint64_t i = 0; i < 2; i++){
        EXPECT_EQ((*bv)[i], ULONG_MAX);
        for(int j = 0; j < 64; j++){
            EXPECT_TRUE(bv->isset(i * 64 + j));
        }
    }
    delete bv;
}

TEST(BitVector, false_constructor){
    BitVector* bv = new BitVector(128, false);
    for(uint64_t i = 0; i < 2; i++){
        EXPECT_EQ((*bv)[i], 0);
        for(uint64_t j = 0; j < 64; j++)
            EXPECT_FALSE(bv->isset(i * 64 + j));
    }
    delete bv;
}


/**
 * Test toggle
 */
TEST(BitVector, toggle_false){
    BitVector* bv = new BitVector(128, true);
    for(uint64_t i = 0; i < 2; i++){
        EXPECT_EQ(bv->operator[](i*64), ULONG_MAX);
        for(uint64_t j = 0; j < 64; j++){
            EXPECT_TRUE(bv->isset(i * 64 + j));
            bv->toggle(i * 64 + j);
        }
    }
    for(uint64_t i = 0; i < 2; i++){
        EXPECT_EQ(bv->operator[](i*64), 0);
        for(uint64_t j = 0; j < 64; j++){
            EXPECT_FALSE(bv->isset(i * 64 + j));
        }
    }
    delete bv;
}

TEST(BitVector, toggle_true){
    BitVector* bv = new BitVector(128, false);
    for(uint64_t i = 0; i < 2; i++){
        EXPECT_EQ(bv->operator[](i*64), 0);
        for(uint64_t j = 0; j < 64; j++){
            EXPECT_FALSE(bv->isset(i * 64 + j));
            bv->toggle(i * 64 + j);
        }
    }
    for(uint64_t i = 0; i < 2; i++){
        EXPECT_EQ(bv->operator[](i*64), ULONG_MAX);
        for(uint64_t j = 0; j < 64; j++){
            EXPECT_TRUE(bv->isset(i * 64 + j));
        }
    }
    delete bv;
}

TEST(BitVector, toggle_half){
    BitVector* bv = new BitVector(128, true);
    for(uint64_t j = 0; j < 64; j++){
        EXPECT_TRUE(bv->isset(j));
        bv->toggle(j);
        EXPECT_FALSE(bv->isset(j));
    }
    EXPECT_EQ(bv->operator[](0), 0);
    EXPECT_EQ(bv->operator[](64), ULONG_MAX);
    delete bv;
}
