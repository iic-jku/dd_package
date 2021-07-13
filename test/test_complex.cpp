/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#include "dd/ComplexNumbers.hpp"

#include "gtest/gtest.h"
#include <memory>

using namespace dd;
using CN = ComplexNumbers;

TEST(DDComplexTest, TrivialTest) {
    auto cn = std::make_unique<ComplexNumbers>();

    unsigned int before_count = cn->cacheCount();

    auto a = cn->getCached(0.5, -0.25);
    auto b = cn->getCached(0.75, 0.75);

    auto r0 = cn->getCached(1.0, -0.125);
    auto r1 = cn->mulCached(a, b);

    cn->lookup(a);
    cn->lookup(b);
    unsigned int between_count = cn->cacheCount();
    // the lookup increases the count in the complex table
    ASSERT_TRUE(before_count < between_count);
    cn->returnToCache(a);
    cn->returnToCache(b);
    cn->returnToCache(r0);
    cn->returnToCache(r1);
    cn->garbageCollect(true);
    // since lookup does not increase the ref count, garbage collection removes the new values
    unsigned int end_count = cn->cacheCount();
    ASSERT_EQ(before_count, end_count);

    EXPECT_NO_THROW(cn->incRef({nullptr, nullptr}));
    EXPECT_NO_THROW(cn->decRef({nullptr, nullptr}));
}

TEST(DDComplexTest, ComplexNumberCreation) {
    auto cn = std::make_unique<ComplexNumbers>();
    ASSERT_EQ(cn->lookup(Complex::zero), Complex::zero);
    ASSERT_EQ(cn->lookup(Complex::one), Complex::one);
    ASSERT_EQ(cn->lookup(1e-14, 0.), Complex::zero);
    ASSERT_EQ(MagEntry::val(cn->lookup(1., 1e-14).mag), 1.);
    ASSERT_EQ(PhaseEntry::val(cn->lookup(1., 1e-14).phase), 0.);
    ASSERT_EQ(MagEntry::val(cn->lookup(1., -1e-14).mag), 1.);
    ASSERT_EQ(PhaseEntry::val(cn->lookup(1., -1e-14).phase), 0.);
    auto c = cn->lookup(1., 1.5);
    std::cout << c << std::endl;
    ASSERT_EQ(MagEntry::val(cn->lookup(c).mag), 1.);
    ASSERT_EQ(PhaseEntry::val(cn->lookup(c).phase), 1.5);
    c = cn->lookup(1., 0.5);
    ASSERT_EQ(MagEntry::val(cn->lookup(c).mag), 1.);
    ASSERT_EQ(PhaseEntry::val(cn->lookup(c).phase), 0.5);
    c = cn->lookup(0.5, -0.5);
    std::cout << c << std::endl;
    ASSERT_EQ(MagEntry::val(cn->lookup(c).mag), 0.5);
    ASSERT_EQ(PhaseEntry::val(cn->lookup(c).phase), 1.5);

    auto e = cn->lookup(1., 2 - 0.25);
    std::cout << e << std::endl;
    std::cout << ComplexValue{1., 0.25} << std::endl;
    cn->magnitudeTable.print();
    cn->magnitudeTable.printStatistics();
    cn->phaseTable.print();
    cn->phaseTable.printStatistics();
}

TEST(DDComplexTest, ComplexNumberArithmetic) {
    auto cn = std::make_unique<ComplexNumbers>();
    auto c  = cn->lookup(1., 0.5);
    auto d  = cn->getCached();
    ComplexNumbers::conj(d, c);
    EXPECT_EQ(MagEntry::val(d.mag), 1.);
    EXPECT_EQ(PhaseEntry ::val(d.phase), 1.5);

    c = cn->lookup(1., 1.75);
    std::cout << "c = " << c.toString(false) << " | phase = " << PhaseEntry::val(c.phase) << " raw phase=" << PhaseEntry::getAlignedPointer(c.phase)->value << "\n";
    ComplexNumbers::neg(d, c);
    std::cout << "d = " << d.toString(false) << " | phase = " << PhaseEntry::val(d.phase) << " raw phase=" << PhaseEntry::getAlignedPointer(d.phase)->value << "\n";
    EXPECT_EQ(MagEntry::val(d.mag), 1.);
    EXPECT_EQ(PhaseEntry::val(d.phase), 0.75);
    c = cn->lookup(dd::SQRT2_2, 0.25);
    ComplexNumbers::incRef(c);
    d = cn->lookup(dd::SQRT2_2, 0.75);
    ComplexNumbers::incRef(d);
    auto e = cn->getTemporary();
    ComplexNumbers::sub(e, c, d);
    ComplexNumbers::decRef(c);
    ComplexNumbers::decRef(d);
    e = cn->lookup(e);
    std::cout << "e.phase=" << PhaseEntry::val(e.phase) << " q=" << +PhaseEntry::getQuadrant(e.phase) << "\n";
    EXPECT_EQ(e, Complex::one);
    auto f = cn->getTemporary();
    ComplexNumbers::div(f, Complex::zero, Complex::one);

    dd::ComplexValue zero{0., 0.};
    dd::ComplexValue one{1., 0.};
    EXPECT_EQ(one + zero, one);
}

TEST(DDComplexTest, NearZeroLookup) {
    auto cn = std::make_unique<ComplexNumbers>();
    auto c  = cn->getTemporary(MagnitudeTable<>::tolerance() / 10., PhaseTable<>::tolerance() / 10.);
    auto d  = cn->lookup(c);
    EXPECT_EQ(d.mag, Complex::zero.mag);
    EXPECT_EQ(d.phase, Complex::zero.phase);
}

TEST(DDComplexTest, SortedBucketsMagnitude) {
    auto     ct  = std::make_unique<MagnitudeTable<>>();
    const fp num = 0.25;

    const std::array<dd::fp, 7> numbers = {
            num + 2. * MagnitudeTable<>::tolerance(),
            num - 2. * MagnitudeTable<>::tolerance(),
            num + 4. * MagnitudeTable<>::tolerance(),
            num,
            num - 4. * MagnitudeTable<>::tolerance(),
            num + 6. * MagnitudeTable<>::tolerance(),
            num + 8. * MagnitudeTable<>::tolerance()};

    const std::size_t the_bucket = ct->hash(num);

    for (auto const& number: numbers) {
        ASSERT_EQ(the_bucket, ct->hash(number));
        ct->lookup(number);
    }

    auto* p = ct->getTable().at(the_bucket);
    ASSERT_NE(p, nullptr);

    dd::fp      last    = std::numeric_limits<dd::fp>::min();
    std::size_t counter = 0;
    while (p != nullptr) {
        ASSERT_LT(last, p->value);
        p = p->next;
        ++counter;
    }
    ct->printStatistics(std::cout);
    EXPECT_EQ(ct->getStatistics().at("lowerNeighbors"), 0);
    EXPECT_EQ(counter, numbers.size());
}

TEST(DDComplexTest, SortedBucketsPhase) {
    auto     ct  = std::make_unique<PhaseTable<>>();
    const fp num = 0.25;

    const std::array<dd::fp, 7> numbers = {
            num + 2. * PhaseTable<>::tolerance(),
            num - 2. * PhaseTable<>::tolerance(),
            num + 4. * PhaseTable<>::tolerance(),
            num,
            num - 4. * PhaseTable<>::tolerance(),
            num + 6. * PhaseTable<>::tolerance(),
            num + 8. * PhaseTable<>::tolerance()};

    const std::size_t the_bucket = ct->hash(num);

    for (auto const& number: numbers) {
        ASSERT_EQ(the_bucket, ct->hash(number));
        ct->lookup(number);
    }

    auto* p = ct->getTable().at(the_bucket);
    ASSERT_NE(p, nullptr);

    dd::fp      last    = std::numeric_limits<dd::fp>::min();
    std::size_t counter = 0;
    while (p != nullptr) {
        ASSERT_LT(last, p->value);
        p = p->next;
        ++counter;
    }
    ct->printStatistics(std::cout);
    EXPECT_EQ(ct->getStatistics().at("lowerNeighbors"), 0);
    EXPECT_EQ(counter, numbers.size());
}

TEST(DDComplexTest, GarbageCollectSomeInBucket) {
    auto cn = std::make_unique<ComplexNumbers>();
    EXPECT_EQ(cn->garbageCollect(), 0);

    const fp num = 0.25;
    cn->lookup(num, 0.0);

    const fp num2 = num + 2. * MagnitudeTable<>::tolerance();
    ComplexNumbers::incRef(cn->lookup(num2, 0.0)); // num2 should be placed in same bucket as num

    auto key  = MagnitudeTable<>::hash(num);
    auto key2 = MagnitudeTable<>::hash(num2);
    ASSERT_EQ(key, key2);

    auto* p = cn->magnitudeTable.getTable()[key];
    EXPECT_NEAR(p->value, num, MagnitudeTable<>::tolerance());

    ASSERT_NE(p->next, nullptr);
    EXPECT_NEAR((p->next)->value, num2, MagnitudeTable<>::tolerance());

    cn->garbageCollect(true); // num should be collected
    EXPECT_NEAR(cn->magnitudeTable.getTable()[key]->value, num2, MagnitudeTable<>::tolerance());
    EXPECT_EQ(cn->magnitudeTable.getTable()[key]->next, nullptr);
}

TEST(DDComplexTest, LookupInNeighbouringBuckets) {
    auto cn = std::make_unique<ComplexNumbers>();

    constexpr std::size_t NBUCKET = MagnitudeTable<>::MASK + 1;

    // lower border of a bucket
    fp bucketBorder = (0.25 * NBUCKET - 0.5) / (NBUCKET - 1);

    // insert a number slightly away from the border
    fp num = bucketBorder + 2 * MagnitudeTable<>::tolerance();
    cn->lookup(num, 0.0);
    auto key = MagnitudeTable<>::hash(num);
    EXPECT_EQ(key, NBUCKET / 4);

    // insert a number barely in the bucket below
    fp num2 = bucketBorder - MagnitudeTable<>::tolerance() / 10;
    cn->lookup(num2, 0.0);
    auto key2 = MagnitudeTable<>::hash(num2);
    EXPECT_EQ(key2, NBUCKET / 4);

    // insert another number in the bucket below a bit farther away from the border
    fp num3 = bucketBorder - 2 * MagnitudeTable<>::tolerance();
    cn->lookup(num3, 0.0);
    auto key3 = MagnitudeTable<>::hash(num3);
    EXPECT_EQ(key3, NBUCKET / 4);

    // insert border number that is too far away from the number in the bucket, but is close enough to a number in the bucket below
    fp   num4 = bucketBorder;
    auto c    = cn->lookup(num4, 0.0);
    auto key4 = MagnitudeTable<>::hash(num4 - MagnitudeTable<>::tolerance());
    EXPECT_EQ(key2, key4);
    EXPECT_NEAR(c.mag->value, num2, MagnitudeTable<>::tolerance());

    // insert a number in the higher bucket
    fp nextBorder = bucketBorder + 1.0 / (NBUCKET - 1);
    fp num5       = nextBorder;
    cn->lookup(num5, 0.0);
    auto key5 = MagnitudeTable<>::hash(num5);
    EXPECT_EQ(key5, NBUCKET / 4 + 1);

    // search for a number in the lower bucket that is ultimately close enough to a number in the upper bucket
    fp   num6 = nextBorder - MagnitudeTable<>::tolerance() / 10;
    auto d    = cn->lookup(num6, 0.0);
    auto key6 = MagnitudeTable<>::hash(num6 + MagnitudeTable<>::tolerance());
    EXPECT_EQ(key5, key6);
    EXPECT_NEAR(d.mag->value, num5, MagnitudeTable<>::tolerance());
}

TEST(DDComplexTest, ComplexValueEquals) {
    ComplexValue a{1.0, 0.0};
    ComplexValue a_tol{1.0 + MagnitudeTable<>::tolerance() / 10, 0.0};
    ComplexValue b{1.0, 0.5};
    EXPECT_TRUE(a.approximatelyEquals(a_tol));
    EXPECT_FALSE(a.approximatelyEquals(b));
}

TEST(DDComplexTest, NumberPrinting) {
    auto cn       = std::make_unique<ComplexNumbers>();
    auto imag     = cn->lookup(1., 0.5);
    auto imag_str = imag.toString(false);
    EXPECT_STREQ(imag_str.c_str(), "1 0.5");
    auto imag_str_formatted = imag.toString(true);
    EXPECT_STREQ(imag_str_formatted.c_str(), "i");

    auto superposition     = cn->lookup(1.0, 0.25);
    auto superposition_str = superposition.toString(false, 3);
    EXPECT_STREQ(superposition_str.c_str(), "1 0.25");
    auto superposition_str_formatted = superposition.toString(true, 3);
    EXPECT_STREQ(superposition_str_formatted.c_str(), "ℯ(iπ/4)");
    auto neg_superposition               = cn->lookup(1.0, -0.25);
    auto neg_superposition_str_formatted = neg_superposition.toString(true, 3);
    EXPECT_STREQ(neg_superposition_str_formatted.c_str(), "ℯ(-iπ/4)");

    std::stringstream ss{};
    ss << std::setprecision(16);
    ComplexValue::printFormatted(ss, dd::SQRT2_2, false);
    EXPECT_STREQ(ss.str().c_str(), "1/√2");
    ss.str("");
    ComplexValue::printFormatted(ss, dd::SQRT2_2, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ/√2)");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.5, false);
    EXPECT_STREQ(ss.str().c_str(), "1/2");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.5, true);
    EXPECT_STREQ(ss.str().c_str(), "i");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.75, false);
    EXPECT_STREQ(ss.str().c_str(), "3/4");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.75, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ 3/4)");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.5 * dd::SQRT2_2, false);
    EXPECT_STREQ(ss.str().c_str(), "1/(2√2)");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.5 * dd::SQRT2_2, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ/(2√2))");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.75 * dd::SQRT2_2, false);
    EXPECT_STREQ(ss.str().c_str(), "3/(4√2)");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.75 * dd::SQRT2_2, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ 3/(4√2))");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.25, false);
    EXPECT_STREQ(ss.str().c_str(), "1/4");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.25, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ/4)");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.25 * dd::SQRT2_2, false);
    EXPECT_STREQ(ss.str().c_str(), "1/(4√2)");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.25 * dd::SQRT2_2, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ/(4√2))");
    ss.str("");

    ComplexValue::printFormatted(ss, dd::PI, false);
    EXPECT_STREQ(ss.str().c_str(), "π");
    ss.str("");
    ComplexValue::printFormatted(ss, dd::PI, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ π)");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.5 * dd::PI, false);
    EXPECT_STREQ(ss.str().c_str(), "π/2");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.5 * dd::PI, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ π/2)");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.25 * dd::PI, false);
    EXPECT_STREQ(ss.str().c_str(), "π/4");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.25 * dd::PI, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ π/4)");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.75 * dd::PI, false);
    EXPECT_STREQ(ss.str().c_str(), "3π/4");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.75 * dd::PI, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ 3π/4)");
    ss.str("");

    ComplexValue::printFormatted(ss, 0.1234567, false);
    EXPECT_STREQ(ss.str().c_str(), "0.1234567");
    ss.str("");
    ComplexValue::printFormatted(ss, 0.1234567, true);
    EXPECT_STREQ(ss.str().c_str(), "ℯ(iπ 0.1234567)");
    ss.str("");
}

TEST(DDComplexTest, MaxRefCountReached) {
    auto cn           = std::make_unique<ComplexNumbers>();
    auto c            = cn->lookup(1., 0.25);
    auto max          = std::numeric_limits<decltype(c.mag->refCount)>::max();
    c.mag->refCount   = max;
    c.phase->refCount = max;
    std::cout.flush();
    std::clog << "Heads up: The following MAXREFCNT warning is part of a passing test.\n";
    CN::incRef(c);
    std::clog.flush();
    EXPECT_EQ(c.mag->refCount, max);
    EXPECT_EQ(c.phase->refCount, max);
}

TEST(DDComplexTest, ComplexTableAllocation) {
    auto        cn     = std::make_unique<ComplexNumbers>();
    std::size_t allocs = cn->magnitudeTable.getAllocations();
    std::cout << allocs << std::endl;
    std::vector<MagnitudeTable<>::Entry*> nums{allocs};
    // get all the numbers that are pre-allocated
    for (auto i = 0U; i < allocs; ++i) {
        nums[i] = cn->magnitudeTable.getEntry();
    }

    // trigger new allocation
    [[maybe_unused]] auto num = cn->magnitudeTable.getEntry();
    EXPECT_EQ(cn->magnitudeTable.getAllocations(), (1. + cn->magnitudeTable.getGrowthFactor()) * allocs);

    // clearing the complex table should reduce the allocated size to the original size
    cn->magnitudeTable.clear();
    EXPECT_EQ(cn->magnitudeTable.getAllocations(), allocs);

    EXPECT_TRUE(cn->magnitudeTable.availableEmpty());
    // obtain entry
    auto entry = cn->magnitudeTable.getEntry();
    // immediately return entry
    cn->magnitudeTable.returnEntry(entry);
    EXPECT_FALSE(cn->magnitudeTable.availableEmpty());
    // obtain the same entry again, but this time from the available stack
    auto entry2 = cn->magnitudeTable.getEntry();
    EXPECT_EQ(entry, entry2);
}

TEST(DDComplexTest, ComplexCacheAllocation) {
    auto cn     = std::make_unique<ComplexNumbers>();
    auto allocs = cn->complexCache.getAllocations();
    std::cout << allocs << std::endl;
    std::vector<Complex> cnums{allocs};
    // get all the cached complex numbers that are pre-allocated
    for (auto i = 0U; i < allocs; i += 1) {
        cnums[i % 2] = cn->getCached();
    }

    EXPECT_EQ(allocs, cn->complexCache.getAllocations());
    // trigger new allocation for obtaining a complex from cache
    [[maybe_unused]] auto cnum = cn->getCached();
    EXPECT_LT(allocs, cn->complexCache.getAllocations());
    EXPECT_EQ(cn->complexCache.getAllocations(), (1 + cn->complexCache.getGrowthFactor()) * allocs);

    // clearing the cache should reduce the allocated size to the original size
    cn->complexCache.clear();
    EXPECT_EQ(cn->complexCache.getAllocations(), allocs);

    // get all the cached complex numbers again
    for (auto i = 0U; i < allocs; i += 1) {
        cnums[i % 2] = cn->getCached();
    }

    // trigger new allocation for obtaining a temporary from cache
    [[maybe_unused]] auto cnumtmp = cn->getTemporary();
    EXPECT_EQ(cn->complexCache.getAllocations(), (1 + cn->complexCache.getGrowthFactor()) * allocs);

    // clearing the unique table should reduce the allocated size to the original size
    cn->complexCache.clear();
    EXPECT_EQ(cn->complexCache.getAllocations(), allocs);
}
