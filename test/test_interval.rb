require 'minitest/autorun'

require 'interval_file_checker'

class TestInterval < MiniTest::Unit::TestCase
  def setup
  end
  
  def test_takes_chr_start_stop
    refute_nil Interval.new(1,1,10)
  end

  def test_b_completely_in_a_overlaps
    a = Interval.new(1,0,10)
    b = Interval.new(1,5,9)
    assert a.overlaps?(b)
    assert b.overlaps?(a)
  end

  def test_overlap_b_starts_in_a
    a = Interval.new(1,0,10)
    b = Interval.new(1,5,15)
    assert a.overlaps?(b)
    assert b.overlaps?(a)
  end

  def test_overlap_small_b
    a = Interval.new(1,1,10)
    b = Interval.new(1,10,10)
    assert a.overlaps?(b)
    assert b.overlaps?(a)
  end

  def test_overlap_different_chr
    a = Interval.new(1,1,10)
    b = Interval.new(2,1,10)
    refute a.overlaps?(b)
    refute b.overlaps?(a)
  end
  
  def test_overlap_close_but_no_cigar
    a = Interval.new(1,1,10)
    b = Interval.new(1,11,20)
    refute a.overlaps?(b)
    refute b.overlaps?(a)
  end

  def test_merge_with_self
    start = 1
    stop = 10
    a = Interval.new(1,start,stop)
    a.merge!(a)
    assert_equal start, a.start
    assert_equal stop, a.stop
  end

  def test_merge_same_chr_a_b
    a = Interval.new(1,1,10)
    b = Interval.new(1,5,15)
    a.merge!(b)
    assert_equal 1, a.start
    assert_equal 15, a.stop
    assert_equal 5, b.start
    assert_equal 15, b.stop
  end

  def test_merge_same_chr_b_a
    a = Interval.new(1,1,10)
    b = Interval.new(1,5,15)
    b.merge!(a)
    assert_equal 1, b.start
    assert_equal 15, b.stop
    assert_equal 1, a.start
    assert_equal 10, a.stop
  end
  
  def test_merge_raises_on_different_chr
    a = Interval.new(1,1,10)
    b = Interval.new(2,1,20)
    assert_raises(ArgumentError) {a.merge!(b)}
  end
end
