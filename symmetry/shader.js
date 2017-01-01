"use strict";

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

(function () {
    let gl;
    let program;
    let vertexshader;
    let fragmentshader;

    let vertexBuffer;
    let showinfo = true;

    let running = false;
    let started = false;
    let error = false;
    let timelast = null;
    let elapsed = 0;
    let frametime = 0;

    let symtype = 0;
    const nsymtypes = 5;
    let ctype = 0;
    const nctypes = 4;
    let ftype = 0;
    const nftypes = 4;

    let P = 4;
    let Q = 4;
    let hplane = 0;
    
    let flags = 0;
    
    let scale = 2.0;
    let rrepeat = 0.0; //1.0/xscale; //1.0/xscale; // One rotation across the screen

    let ulimit = 1.0;
    let vlimit = 1.0;
    let uscale = 1.0;
    let uxfact = 0.0;
    let uyfact = 0.0;
    let utfact = 0.0;

    let vscale = 1.0;
    let vxfact = 0.0;
    let vyfact = 0.0;
    let vtfact = 0.0;

    let xoffset = 0;
    let yoffset = 0;
    let uoffset = 0;
    let voffset = 0;
    let time = 0;
    
    // Initialize WebGL, returning the GL context or null if
    // WebGL isn't available or could not be initialized.
    function initWebGL(canvas,attributes) {
        let gl = null;
        try {
            // Could have "webgl-experimental" here instead
            gl = canvas.getContext("experimental-webgl")
        }
        catch(e) {
            console.log(e)
        }

        // If we don't have a GL context, give up now
        if (!gl) alert("Unable to initialize WebGL");
        error = true;
        return gl;
    }

    // Initialize the shaders, so WebGL knows how to light our scene.
    function initShaders(vshader,fshader) {
        let vertexShader = makeShader(vshader.source,gl.VERTEX_SHADER);
        let fragmentShader = makeShader(fshader.source,gl.FRAGMENT_SHADER);
        if (!vertexShader || !fragmentShader) return false;
        let program = gl.createProgram();
        gl.attachShader(program, vertexShader);
        gl.attachShader(program, fragmentShader);
        gl.linkProgram(program);
        gl.validateProgram(program); // Check all well
        // If creating the shader program failed, alert
        if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
            alert("Unable to initialize the shader program: " +
                  gl.getProgramInfoLog(program));
            error = true;
            return false;
        }
        return program;
    }

    function makeShader(source, shadertype) {
        const shader = gl.createShader(shadertype);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            alert("Error for " + shadertype + ": " + gl.getShaderInfoLog(shader));
            error = true;
            return null;
        }
        return shader;
    }

    // Loads a shader program by scouring the current document,
    // looking for a script with the specified ID.
    function getShader(id) {
        const shaderelement = document.getElementById(id);
        // Didn't find an element with the specified ID; abort.
        if (!shaderelement) {
            return null;
        }

        // Walk through the source element's children, building the
        // shader source string.
        let source = "";
        let currentChild = shaderelement.firstChild;

        while(currentChild) {
            if (currentChild.nodeType == 3) {
                source += currentChild.textContent;
            }
            currentChild = currentChild.nextSibling;
        }

        // Now figure out what type of shader script we have,
        // based on its MIME type.
        let shadertype;
        if (shaderelement.type == "x-shader/x-fragment") {
            shadertype = gl.FRAGMENT_SHADER;
        } else if (shaderelement.type == "x-shader/x-vertex") {
            shadertype = gl.VERTEX_SHADER;
        } else {
            return null;  // Unknown shader type
        }
        return makeShader(source,shadertype);
    }

    let resourcesLoading = 0;
    function resourceLoaded() {
        resourcesLoading--;
        if (resourcesLoading == 0) {
            program = initShaders(vertexshader,fragmentshader);
            if (program) {
                initBuffers();
                window.addEventListener("keypress",keypressHandler,false);
                window.addEventListener("resize",function () {
                    if (started && !running) requestAnimationFrame(drawScene);
                });
                // Use window param explicitly here.
                window.requestAnimationFrame(drawScene);
                started = true;
            }
        }
    }
    
    function initTexture(filename,format,unit) {
        let texture = gl.createTexture();
        // new Image() makes an HTMLImageElement()
        // We store this as a property of the texture object.
        texture.image = new Image();
        texture.image.onload = function() {
            handleLoadedTexture(texture,format,unit)
        }
        texture.image.src = filename;
        resourcesLoading++;
        return true;
    }

    function handleLoadedTexture(texture,format,unit) {
        // Bind the texture with the sampler unit.
        // Make the texture active (on unit 3 - as per uniform setting above)
        gl.activeTexture(unit);
        // "Subsequent calls to bindTexture will bind the
        // texture to the currently active unit"
        gl.bindTexture(gl.TEXTURE_2D, texture);
        gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
        //gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texture.image);
        //gl.texImage2D(gl.TEXTURE_2D, 0, gl.LUMINANCE_ALPHA, gl.LUMINANCE_ALPHA, gl.UNSIGNED_BYTE, texture.image);
        gl.texImage2D(gl.TEXTURE_2D, 0, format, format, gl.UNSIGNED_BYTE, texture.image);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        //gl.hint(gl.GENERATE_MIPMAP_HINT, gl.NICEST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
        // For patterns that line up at the edge, may not want mirroring.
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.MIRRORED_REPEAT);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.MIRRORED_REPEAT);
        gl.generateMipmap(gl.TEXTURE_2D);
        // Only once the textures are loaded do we try and draw
        resourceLoaded();
    }

    function initshader(url) {
        const request = new XMLHttpRequest();
        const shader = {}
        request.open("GET", url);
        request.onreadystatechange = function() {
            //console.log("onreadystatechange", request.readyState);
            if (request.readyState === 4) {
                //console.log(request.status, request.getResponseHeader("Content-Type"));
                shader.source = request.responseText;
                resourceLoaded();
            }
        }
        resourcesLoading++;
        request.send(null); // No body
        return shader;
    }

    // Allocate and initialize storage on the GPU for the scene data.
    // Later these buffers will be associated with per-vertex attributes.
    function initBuffers() {
        const vertices = [
            1.0, 1.0, 0.0,
           -1.0, 1.0, 0.0,
            1.0,-1.0, 0.0,
           -1.0,-1.0, 0.0
        ];
        const colors = [
            1.0,  1.0,  1.0,  1.0,    // white
            1.0,  0.0,  0.0,  1.0,    // red
            0.0,  1.0,  0.0,  1.0,    // green
            0.0,  0.0,  1.0,  1.0     // blue
        ];
        // Texture coordinates
        const uvs = [
            1.0, 1.0,
            0.0, 1.0,
            1.0, 0.0,
            0.0, 0.0
        ];
        const vertexArray = new Float32Array(vertices.length + colors.length + uvs.length);
        const VSIZE = 3, CSIZE = 4, UVSIZE = 2;
        let off = 0;
        for (let i = 0; i < vertices.length/VSIZE; i++) {
            for (let j = 0; j < VSIZE; j++) {
                vertexArray[off++] = vertices[VSIZE*i+j];
            }
            for (let j = 0; j < CSIZE; j++) {
                vertexArray[off++] = colors[CSIZE*i+j];
            }
            for (let j = 0; j < UVSIZE; j++) {
                vertexArray[off++] = uvs[UVSIZE*i+j];
            }
        }
        console.assert(off == vertexArray.length);
        
        vertexBuffer = gl.createBuffer();
        // Select as ARRAY_BUFFER target
        gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
        // Now pass the list of vertices into WebGL to build the shape. We
        // do this by creating a Float32Array from the JavaScript array,
        // then use it to fill the current vertex buffer.

        // If we wanted to modify the vertex data, we would
        // modify vertexArray and call gl.bufferData again (or
        // use gl.bufferSubData(target,offset,data)) and we
        // would want gl.DYNAMIC_DRAW as usage as well.
        gl.bufferData(gl.ARRAY_BUFFER, vertexArray, gl.STATIC_DRAW);
        gl.bindBuffer(gl.ARRAY_BUFFER, null);
        return true;
    }

    const helpstring =
          "&lt;space&gt;: scroll s/S: scale r/R: rrepeat t/T: utfact " +
          "u/U: uxfact v/V: vytype x/X: symtype c/C: ctype f/F: ftype 1: hexagonal " +
          "?: info !: mkurl"

    function keypressHandler(event) {
        if (!event.ctrlKey) {
            // Ignore event if control key pressed.
            var c = String.fromCharCode(event.charCode)
            var handled = true;
            switch(c) {
            case ' ':
                running = !running;
                if (started && running) {
                    // If we are now running, start animating.
                    requestAnimationFrame(drawScene);
                    timelast = null;
                }
                break;
            case '1': flags ^= 1; break;
            case '2': flags ^= 2; break;
            case '3': flags ^= 4; break;
            case '4': flags ^= 8; break;
            case '5': flags ^= 16; break;
            //case '6': flags ^= 32; break;
            case 'p': P = (P+1)%16; break;
            case 'P': P = (P+15)%16; break;
            case 'q': Q = (Q+1)%16; break;
            case 'Q': Q = (Q+15)%16; break;
            case 's': scale *= 1.1; break;
            case 'S': scale /= 1.1; break;
            case 'r': rrepeat += 0.05; break;
            case 'R': rrepeat -= 0.05; break;
            case 't': utfact += 0.05; break;
            case 'T': utfact -= 0.05; break;
            case 'u': uxfact += 0.1; break;
            case 'U': uxfact -= 0.1; break;
            case 'v': vyfact += 0.1; break;
            case 'V': vyfact -= 0.1; break;
                // How this works is all wrong!
            case 'x': xoffset += 0.1; break;
            case 'X': xoffset -= 0.1; break;
            case 'y': yoffset += 0.1; break;
            case 'Y': yoffset -= 0.1; break;
            case 'a': symtype = (symtype+1)%nsymtypes; break;
            case 'A': symtype = (symtype+nsymtypes-1)%nsymtypes; break;
            case 'c': ctype = (ctype+1)%nctypes; break;
            case 'C': ctype = (ctype+nctypes-1)%nctypes; break;
            case 'f': ftype = (ftype+1)%nftypes; break;
            case 'F': ftype = (ftype+nftypes-1)%nftypes; break;
            case 'h': hplane = (hplane+1)%5; break;
            case 'H': hplane = (hplane+4)%5; break;
            case '?': showinfo = !showinfo; setinfo(); break;
            case '!': alert(mkurl()); break;
            default:
                handled = false;
                break;
            }
            if (started && !running) requestAnimationFrame(drawScene);
            if (handled) event.preventDefault();
        }
    }

    function makeflags() {
        return flags + (P<<5) + (Q<<(5+4)) + (hplane<<(5+4+4));
    }
    function initProgram(delta) {
        gl.useProgram(program);
        
        let xscale = scale;
        let yscale = scale*(gl.canvas.height/gl.canvas.width);
        let xscroll = 0.0; //scale*0.05; // Scrolling speed scale*1.0 is whole width in 1 second
        let yscroll = 0.0;

        // Finally associate the buffer with the vertexPositionAttribute
        // as defined in the shader. Is there any reason to delay doing this?
        // Probably only valid to do when the correct program is active.
        // We could autocalculate the stride and offset.
        const buffers = [ { name: "aVertexPosition",
                            buffer: vertexBuffer,
                            size: 3,
                            type: gl.FLOAT,
                            stride: 9*4, // bytes
                            offset: 0*4  // bytes
                          },
                          { name: "aVertexColor",
                            buffer: vertexBuffer,
                            size: 4,
                            type: gl.FLOAT,
                            stride: 9*4,
                            offset: 3*4
                          },
                          { name: "aVertexUV",
                            buffer: vertexBuffer,
                            size: 2,
                            type: gl.FLOAT,
                            stride: 9*4,
                            offset: 7*4
                          } ];
        for (let i = 0; i < buffers.length; i++) {
            let b = buffers[i]
            const index = gl.getAttribLocation(program, b.name);
            if (index >= 0) {
                gl.bindBuffer(gl.ARRAY_BUFFER, b.buffer);
                gl.enableVertexAttribArray(index);
                gl.vertexAttribPointer(index, b.size, b.type, false, b.stride, b.offset);
            }
        }
        // locations are fixed for a given program
        const uniforms = [
            { location: gl.getUniformLocation(program, "uSampler"), value: 1 },
            { location: gl.getUniformLocation(program, "uNoise"), value: 2 },
        ];

        // potentially change the values each time around
        for (let i = 0; i < uniforms.length; i++) {
            let u = uniforms[i];
            gl.uniform1i(u.location, u.value);
        }

        // These get turned in to 32 bit floats on the GPU so
        // restrict the range to avoid loss of precision.
        // xoffset = (xoffset + xscroll*delta) % 100;
        // yoffset = 0; // (yoffset + yscroll*delta) % 100;
        
        uoffset = (uoffset + utfact*delta) % 100;
        voffset = (voffset + vtfact*delta) % 100;

        time = (time + delta) % 100;;
        
        gl.uniform4f(gl.getUniformLocation(program,"params1"),
                     xscale, yscale, xoffset, yoffset);
        gl.uniform4f(gl.getUniformLocation(program,"params2"),
                     ulimit,vlimit,rrepeat,0);

        let sUniform = gl.getUniformLocation(program, "uVScale")
        gl.uniform1f(sUniform, gl.canvas.height/gl.canvas.width)
        const ufactUniform = gl.getUniformLocation(program, "ufact");
        gl.uniform4f(ufactUniform,uscale,uxfact,uyfact,uoffset);
        const vfactUniform = gl.getUniformLocation(program, "vfact");
        gl.uniform4f(vfactUniform,vscale,vxfact,vyfact,voffset);

        const ABUniform = gl.getUniformLocation(program, "AB");
        const CDUniform = gl.getUniformLocation(program, "CD");
        gl.uniform4f(ABUniform,Math.cos(0.1*time),Math.sin(0.11*time),1,0);
        gl.uniform4f(CDUniform,1,0,0,0);

        gl.uniform1i(gl.getUniformLocation(program, "uFlags"), makeflags());
        
        const aUniform = gl.getUniformLocation(program, "a");
        switch(ctype) {
        case 0:
            gl.uniform4f(aUniform,0.25,0.25,0.25,0.25);
            break;
        case 1:
            gl.uniform4f(aUniform,0.4,0.4,0.1,0.1);
            break;
        case 2:
            gl.uniform4f(aUniform,1.0,0.0,0.0,0.0);
            break;
        case 3:
            gl.uniform4f(aUniform,
                         0.4*Math.sin(0.1*time), 0.4*Math.cos(0.1*time),
                         0.1*Math.sin(0.11*time),0.1*Math.cos(0.11*time));
            break;
        default:
            console.assert(0);
        }

        const mn0Uniform = gl.getUniformLocation(program, "mn0");
        const mn1Uniform = gl.getUniformLocation(program, "mn1");
        switch(symtype) {
            // Should use an array for this.
        case 0:
            gl.uniform4i(mn0Uniform,0,1,1,0); gl.uniform4i(mn1Uniform,1,2,-2,1); // Generic rotation
            break;
        case 1:
            gl.uniform4i(mn0Uniform,0,1,1,0); gl.uniform4i(mn1Uniform,1,2,2,1); // p31m
            break;
        case 2:
            gl.uniform4i(mn0Uniform,0,1,-1,0); gl.uniform4i(mn1Uniform,1,2,-2,-1); // p3m1
            break;
        case 3:
            gl.uniform4i(mn0Uniform,0,1,0,-1); gl.uniform4i(mn1Uniform,1,2,-1,-2); // p6
            break;
        case 4:
            gl.uniform4i(mn0Uniform,0,1,0,-1); gl.uniform4i(mn1Uniform,1,0,-1,0); // p6m
            break;
        default:
            console.log(symtype);
            console.assert(0);
        }
    }

    function renderScene(delta) {
        initProgram(delta);
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    }

    // Draw the scene.
    function drawScene(timenow) {
        const canvas = gl.canvas
        let delta = 0;
        console.assert(started);
        if (running) {
            //setTimeout(function () { requestAnimationFrame(drawScene) }, 30);
            requestAnimationFrame(drawScene);
            if (timelast == null) {
                timelast = timenow;
            } else {
                frametime = 0.99*frametime+0.01*(timenow-timelast)
            }
            delta = timenow-timelast;
            timelast = timenow;
        }
        if (showinfo) {
            let s = "";
            s += "Framerate: " + ((1000/frametime)|0);
            s += " symtype: " + symtype;
            s += " ctype: " + ctype;
            s += " ftype: " + ftype;
            s += " scale: " + scale.toFixed(2);
            s += " rrepeat: " + rrepeat.toFixed(2);
            s += " uxfact: " + uxfact.toFixed(2);
            s += " vyfact: " + vyfact.toFixed(2);
            s += " utfact: " + utfact.toFixed(2);
            info.innerHTML = s;
        }
        
        // Check for resize
        //let width = gl.drawingBufferWidth;
        //let height = gl.drawingBufferHeight;
        let width = canvas.clientWidth;
        let height = canvas.clientHeight;
        if (width != canvas.width || height != canvas.height) {
            canvas.width = width;
            canvas.height = height;
        }
        gl.viewport(0,0,canvas.width,canvas.height);
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        renderScene(delta/1000); // Time since last render in seconds
    }

    var imgname = "color.jpg";
    function mkurl() {
        let s = "?";
        s += "img=" + imgname;
        if (hexagonal) s += "&hexagonal";
        else s += "&square";
        if (scale != 2.0) s += "&scale=" + scale.toFixed(2);
        if (symtype != 0) s += "&symtype=" + symtype;
        if (ctype != 0) s += "&ctype=" + ctype;
        if (ftype != 0) s += "&ftype=" + ftype;
        if (rrepeat != 0) s += "&rrepeat=" + rrepeat.toFixed(2);
        if (uxfact != 0) s += "&uxfact=" + uxfact.toFixed(2);
        if (uyfact != 0) s += "&uyfact=" + uyfact.toFixed(2);
        if (utfact != 0) s += "&utfact=" + utfact.toFixed(2);
        if (uscale != 1.0) s += "&uscale=" + uscale.toFixed(2);
        if (ulimit != 1.0) s += "&ulimit=" + ulimit.toFixed(2);

        if (vxfact != 0) s += "&vxfact=" + vxfact.toFixed(2);
        if (vyfact != 0) s += "&vyfact=" + vyfact.toFixed(2);
        if (vtfact != 0) s += "&vtfact=" + vtfact.toFixed(2);
        if (vscale != 1.0) s += "&vscale=" + vscale.toFixed(2);
        if (vlimit != 1.0) s += "&vlimit=" + vlimit.toFixed(2);
        return s;
    }
    function setinfo() {
        if (showinfo) {
            info.style.display = 'block';
            help.style.display = 'block';
        } else {
            info.style.display = 'none';
            help.style.display = 'none';
        }
    }
        
    window.runoncanvas = function(canvas,vsfile,fsfile) {
        help.innerHTML = helpstring;
        setinfo();
        var options = window.location.search;
        if (options.length > 0) {
            // Strip off leading '?'
            options = options.slice(1).split('&');
            options.forEach(function(arg) {
                let matches;
                if (matches = arg.match(/^img=(.+)$/)) {
                    imgname = matches[1];
                } else if (matches = arg.match(/^hexagonal$/)) {
                    flags |= 1;
                } else if (matches = arg.match(/^square$/)) {
                    flags &= ~1;
                } else if (matches = arg.match(/^scale=([\d.-]+)$/)) {
                    scale = Number(matches[1]);
                } else if (matches = arg.match(/^symtype=([\d]+)$/)) {
                    symtype = Number(matches[1])%nsymtypes;
                } else if (matches = arg.match(/^ctype=([\d]+)$/)) {
                    ctype = Number(matches[1])%nctypes;
                } else if (matches = arg.match(/^ftype=([\d]+)$/)) {
                    ftype = Number(matches[1])%nftypes;
                } else if (matches = arg.match(/^rrepeat=([\d.-]+)$/)) {
                    rrepeat = Number(matches[1]);
                } else if (matches = arg.match(/^uxfact=([\d.-]+)$/)) {
                    uxfact = Number(matches[1]);
                } else if (matches = arg.match(/^uyfact=([\d.-]+)$/)) {
                    uyfact = Number(matches[1]);
                } else if (matches = arg.match(/^vxfact=([\d.-]+)$/)) {
                    vxfact = Number(matches[1]);
                } else if (matches = arg.match(/^vyfact=([\d.-]+)$/)) {
                    vyfact = Number(matches[1]);
                } else if (matches = arg.match(/^utfact=([\d.-]+)$/)) {
                    utfact = Number(matches[1]);
                } else if (matches = arg.match(/^vtfact=([\d.-]+)$/)) {
                    vtfact = Number(matches[1]);
                } else if (matches = arg.match(/^uscale=([\d.-]+)$/)) {
                    uscale = Number(matches[1]);
                } else if (matches = arg.match(/^vscale=([\d.-]+)$/)) {
                    vscale = Number(matches[1]);
                } else if (matches = arg.match(/^ulimit=([\d.-]+)$/)) {
                    ulimit = Number(matches[1]);
                } else if (matches = arg.match(/^vlimit=([\d.-]+)$/)) {
                    vlimit = Number(matches[1]);
                } else {
                    console.log("Ignoring parameter '" + arg + "'");
                }
            });
        }

        gl = initWebGL(canvas);      // Initialize the GL context
        // Only continue if WebGL is available and working
        if (gl) {
            gl.clearColor(0.0, 0.0, 0.0, 1.0);  // Set clear color to black and fully opaque
            initTexture("../images/" + imgname, gl.RGBA, gl.TEXTURE1);
            vertexshader = initshader(vsfile);
            fragmentshader = initshader(fsfile);
            setTimeout(function(){
                if (!started && !error) {
                    alert("Page load timed out: " + resourcesLoading);
                }
            }, 5000);
        }
    }
})();
