"use strict"

// The MIT License (MIT)

// Copyright (c) 2017 Matthew Arcus

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Shouldn't be at top level
function makeiframes(urls) {
    for (let i = 0; i < urls.length; i++) {
        const url = urls[i];
        const iframe = document.createElement("iframe");
        iframe.src = url;
        iframe.width = "100%";
        iframe.height = "600";
        const para = document.createElement("p");
        const link = document.createElement("a");
        para.appendChild(link);
        link.href = url;
        const text = document.createTextNode(url);
        link.appendChild(text);
        document.body.appendChild(iframe);
        document.body.appendChild(para);
    }
}

function hsva(h,s,v,a) {
    // 1,3,5 are pure R,G,B
    h %= 1;
    if (a == undefined) a = 1
    h *= 6;
    const f = h%1
    let r = 0, g = 0, b = 0
    if (h < 1) {
        r = 1; b = 1-f;
    } else if (h < 2) {
        r = 1; g = f;
    } else if (h < 3) {
        g = 1; r = 1-f;
    } else if (h < 4) {
        g = 1; b = f;
    } else if (h < 5) {
        b = 1; g = 1-f;
    } else {
        b = 1; r = f;
    }
    r = v*(s*r+1-s);
    g = v*(s*g+1-s);
    b = v*(s*b+1-s);
    return [r,g,b,a]
}

function test() {
    let N = 10;
    for (let i = 0; i <= N; i++) {
        for (let j = 0; j <= N; j++) {
            for (let k = 0; k <= N; k++) {
                console.log(hsva(i/N,j/N,k/N,1));
            }
        }
    }
}
//test()
