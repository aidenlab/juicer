#!/usr/bin/awk -f
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# Deduping script, checks that positions are within wobble or not
# Operates on SAM files with rt and cb flags set
# Juicer version 2.0
#
# Streaming dedup: O(n*k) time, O(k) memory
# Streams: when a read enters a group, immediately check it against
# the kept-window and print it (as dup or not). No group buffer needed.
# Memory = O(k) where k = number of unique positions kept (tiny).
#
# For small groups (<=100), still buffer+O(n^2) for exact parity with original.
# For large groups, stream with kept-window.

function markduplicate(input,    _p1,_p2,_flag) {
    _p1 = index(input, "\t");
    _p2 = index(substr(input, _p1+1), "\t");
    _flag = substr(input, _p1+1, _p2-1) + 0;
    if (int(_flag / 1024) % 2 == 0)
        print substr(input, 1, _p1) (_flag + 1024) substr(input, _p1+_p2);
    else
        print input;
}

# Print current read's lines (not a dup)
function print_cur() {
    print cur_l1;
    if (cur_n >= 2) print cur_l2;
    if (cur_n >= 3) print cur_l3;
    if (cur_n >= 4) print cur_l4;
}
# Mark current read's lines as duplicate and print
function mark_cur() {
    markduplicate(cur_l1);
    if (cur_n >= 2) markduplicate(cur_l2);
    if (cur_n >= 3) markduplicate(cur_l3);
    if (cur_n >= 4) markduplicate(cur_l4);
    count_dups++;
}

# Check current read against kept window. Return 1 if dup, 0 if not.
function check_dup(p1val, p2val,    _k, _d1, _d2) {
    for (_k = nkept - 1; _k >= wstart; _k--) {
        _d1 = p1val - kept_p1[_k]; if (_d1 < 0) _d1 = -_d1;
        if (_d1 > wobble1) break;
        _d2 = p2val - kept_p2[_k]; if (_d2 < 0) _d2 = -_d2;
        if ((_d1 <= wobble1 && _d2 <= wobble2) || (_d2 <= wobble1 && _d1 <= wobble2))
            return 1;
    }
    return 0;
}

# Add a read to the group in STREAMING mode: check, print, update kept window
function stream_add(p1val, p2val) {
    # Advance window start past reads whose pos1 is too far
    while (wstart < nkept && kept_p1[wstart] < p1val - wobble1) wstart++;

    if (check_dup(p1val, p2val)) {
        mark_cur();
    } else {
        print_cur();
        kept_p1[nkept] = p1val;
        kept_p2[nkept] = p2val;
        nkept++;
    }
}

# Start a new streaming group with the current read as first (non-dup) entry
function stream_start(p1val, p2val) {
    # Don't need to delete kept arrays — overwritten by index, wstart handles window
    nkept = 0; wstart = 0;
    print_cur();
    kept_p1[0] = p1val;
    kept_p2[0] = p2val;
    nkept = 1;
}

# === Small group buffered approach (unchanged from original) ===
function buf_add(p1val, p2val) {
    grp_l1[gi]=cur_l1; grp_n[gi]=cur_n;
    if (cur_n>=2) grp_l2[gi]=cur_l2;
    if (cur_n>=3) grp_l3[gi]=cur_l3;
    if (cur_n>=4) grp_l4[gi]=cur_l4;
    gpos1[gi]=p1val; gpos2[gi]=p2val; gi++;
}
function buf_flush() {
    if (gi == 1) {
        print grp_l1[0];
        if (grp_n[0] >= 2) print grp_l2[0];
        if (grp_n[0] >= 3) print grp_l3[0];
        if (grp_n[0] >= 4) print grp_l4[0];
    } else if (gi > 1) {
        for (j=0; j<gi; j++) {
            if (!(j in dups)) {
                for (k=j+1; k<gi; k++) {
                    _d1=gpos1[k]-gpos1[j]; if(_d1<0)_d1=-_d1;
                    _d2=gpos2[k]-gpos2[j]; if(_d2<0)_d2=-_d2;
                    if ((_d1<=wobble1 && _d2<=wobble2) || (_d2<=wobble1 && _d1<=wobble2)) dups[k]++;
                    if (_d1>wobble1) break;
                }
            }
        }
        for (j=0; j<gi; j++) {
            if (j in dups) {
                count_dups++;
                markduplicate(grp_l1[j]);
                if (grp_n[j] >= 2) markduplicate(grp_l2[j]);
                if (grp_n[j] >= 3) markduplicate(grp_l3[j]);
                if (grp_n[j] >= 4) markduplicate(grp_l4[j]);
            } else {
                print grp_l1[j];
                if (grp_n[j] >= 2) print grp_l2[j];
                if (grp_n[j] >= 3) print grp_l3[j];
                if (grp_n[j] >= 4) print grp_l4[j];
            }
        }
        delete dups;
    }
}
# Transition from buffered to streaming mode when group exceeds threshold
function buf_to_stream() {
    # Re-process buffered entries through the kept-window
    nkept = 0; wstart = 0;
    for (j=0; j<gi; j++) {
        # Check each buffered entry against kept window
        _is_dup = 0;
        for (_k = nkept - 1; _k >= wstart; _k--) {
            _d1 = gpos1[j] - kept_p1[_k]; if (_d1 < 0) _d1 = -_d1;
            if (_d1 > wobble1) break;
            _d2 = gpos2[j] - kept_p2[_k]; if (_d2 < 0) _d2 = -_d2;
            if ((_d1<=wobble1 && _d2<=wobble2) || (_d2<=wobble1 && _d1<=wobble2)) { _is_dup=1; break }
        }
        if (_is_dup) {
            count_dups++;
            markduplicate(grp_l1[j]);
            if (grp_n[j] >= 2) markduplicate(grp_l2[j]);
            if (grp_n[j] >= 3) markduplicate(grp_l3[j]);
            if (grp_n[j] >= 4) markduplicate(grp_l4[j]);
        } else {
            print grp_l1[j];
            if (grp_n[j] >= 2) print grp_l2[j];
            if (grp_n[j] >= 3) print grp_l3[j];
            if (grp_n[j] >= 4) print grp_l4[j];
            kept_p1[nkept] = gpos1[j];
            kept_p2[nkept] = gpos2[j];
            nkept++;
        }
    }
    streaming = 1;
}

BEGIN {
    gi = 0; cur_n = 0; pname = ""; _started = 0;
    streaming = 0;  # 0 = buffered mode, 1 = streaming mode
    STREAM_THRESHOLD = 100;
    if (length(nowobble)>0) { wobble1=0; wobble2=0 }
    else if (length(wobble1)==0 && length(wobble2)==0) { wobble1=4; wobble2=4 }
    neg_wobble1 = -wobble1;
    OFS = "\t";
}
{
    if (NF < 11) { print; next }
    _cbtag = $(NF-1);
    if (substr(_cbtag, 1, 4) != "cb:Z") { _cbtag = $NF; if (substr(_cbtag, 1, 4) != "cb:Z") { print; next } }

    if (pname == $1) {
        cur_n++;
        if (cur_n == 2) cur_l2 = $0;
        else if (cur_n == 3) cur_l3 = $0;
        else cur_l4 = $0;
        next;
    }

    # NEW read name — process previous read
    if (_started) {
        _d = cc7 - p7;
        if (_d >= neg_wobble1 && _d <= wobble1 && cc1==p1 && cc2==p2 && cc3==p3 && cc4==p4 && cc5==p5 && cc6==p6) {
            # Previous read belongs in current group
            if (streaming) {
                stream_add(cc7, cc8);
            } else {
                buf_add(cc7, cc8);
                if (gi >= STREAM_THRESHOLD) {
                    buf_to_stream();
                }
            }
        } else {
            # Previous read starts a new group — flush current group
            if (streaming) {
                # Streaming group is already flushed (printed as we went)
                streaming = 0;
            } else {
                buf_flush();
            }
            # Start new group with previous read
            gi = 0;
            buf_add(cc7, cc8);
        }
        p1=cc1; p2=cc2; p3=cc3; p4=cc4; p5=cc5; p6=cc6; p7=cc7; p8=cc8;
    }

    # Start accumulating new read
    pname = $1;
    cur_l1 = $0;
    cur_n = 1;
    split(substr(_cbtag, 6), _cba, "_");
    cc1=_cba[1]; cc2=_cba[2]; cc3=_cba[3]+0; cc4=_cba[4]+0;
    cc5=_cba[5]; cc6=_cba[6]; cc7=_cba[7]+0; cc8=_cba[8]+0;
    _started = 1;
}
END {
    if (_started) {
        _d = cc7 - p7;
        if (_d >= neg_wobble1 && _d <= wobble1 && cc1==p1 && cc2==p2 && cc3==p3 && cc4==p4 && cc5==p5 && cc6==p6) {
            if (streaming) {
                stream_add(cc7, cc8);
            } else {
                buf_add(cc7, cc8);
                if (gi >= STREAM_THRESHOLD) buf_to_stream();
            }
        } else {
            if (streaming) { streaming=0 } else { buf_flush() }
            gi = 0;
            buf_add(cc7, cc8);
        }
        if (streaming) {
            # done, already printed
        } else {
            buf_flush();
        }
    }
}

