"use strict";
/* TODO
Use fetch instead of XMLHttpRequest
Periodic retry of shader downloads (or just on demand)
Multiple textures (make configurable)
Shadertoy style keyboard texture (or uniform)
Popup menu for config
 */
(function () {
    let gl;
    let program;
    let vertexshader;
    let fragmentshader;

    let buffers;
    let vertexBuffer;
    let showinfo = false;

    let running = false;
    let started = false;
    let error = false;
    let timelast = null;
    let frametime = 0;

    let itime = 0;
    let progressive = false;

    let mousex = 0;
    let mousey = 0;
    
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
        if (!gl.getProgramParameter(program, gl.LINK_STATUS) || gl.isContextLost()) {
            alert("Unable to initialize the shader program: " +
                  gl.getProgramInfoLog(program));
            error = true;
            return false;
        }
        return program;
    }

    function makeShader(source, shadertype) {
        const shader = gl.createShader(shadertype);
        if (!source.match("^#version")) {
            source = "#version 300 es\n#define HEADER\n#line 1\n" + source;
        }
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            alert("Error for " + shadertype + ": " + gl.getShaderInfoLog(shader));
            error = true;
            return null;
        }
        return shader;
    }

    function startreloadtimer(shader) {
        setTimeout(function () {
            console.log("Reloadtimer: " + shader.url + " " + shader.lastmodified);
            startreloadtimer(shader);
        },
                   1000);
    }
    
    let framenumber = 0;
    let resourcesLoading = 0;
    let cubeTexture = null;
    function resourceLoaded() {
        resourcesLoading--;
        if (resourcesLoading == 0 && !error) {
            //startreloadtimer(fragmentshader);
            //console.log("Resources loaded");
            program = initShaders(vertexshader,fragmentshader);
            if (program) {
                if (cubeTexture) {
                    gl.activeTexture(gl.TEXTURE3);
                    gl.bindTexture(gl.TEXTURE_CUBE_MAP, cubeTexture);
                    //gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
                    gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
                    gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
                    gl.hint(gl.GENERATE_MIPMAP_HINT, gl.NICEST);
                    gl.generateMipmap(gl.TEXTURE_CUBE_MAP);
                }
                initBuffers();

	        window.addEventListener( 'mousedown', onMouseDown, false );
	        window.addEventListener( 'wheel', onMouseWheel, false );

	        window.addEventListener( 'touchstart', onTouchStart, false );
	        window.addEventListener( 'touchmove', onTouchMove, false );
	        window.addEventListener( 'touchend', onTouchEnd, false );

                window.addEventListener("keydown",keydownHandler,false);
                window.addEventListener("keyup",keyupHandler,false);
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
    
    function initCubeTexture(filename,format,type) {
        // new Image() makes an HTMLImageElement()
        // We store this as a property of the texture object.
        let image = new Image();
        image.onload = function() {
            handleLoadedCubeTexture(image,format,type,filename);
        }
        image.src = filename;
        resourcesLoading++;
        return true;
    }

    function handleLoadedCubeTexture(image,format,type,fname) {
        gl.bindTexture(gl.TEXTURE_CUBE_MAP, cubeTexture);
        gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
        gl.texImage2D(type, 0, format, format, gl.UNSIGNED_BYTE, image);
        //gl.texParameteri(gl.TEXTURE_CUBE_MAP, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
        resourceLoaded();
    }

    function initTexture(filename,format,unit) {
        console.log("Init texture: " + filename);
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
        console.log("Texture loaded: " + texture + " " + unit);
        gl.activeTexture(unit); // Activate unit
        // "Subsequent calls to bindTexture will bind the texture to the currently active unit"
        gl.bindTexture(gl.TEXTURE_2D, texture);
        gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
        //gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texture.image);
        //gl.texImage2D(gl.TEXTURE_2D, 0, gl.LUMINANCE_ALPHA, gl.LUMINANCE_ALPHA, gl.UNSIGNED_BYTE, texture.image);
        gl.texImage2D(gl.TEXTURE_2D, 0, format, format, gl.UNSIGNED_BYTE, texture.image);
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST_MIPMAP_LINEAR); // This is default anyway
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
        // For patterns that line up at the edge, may not want mirroring.
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.MIRRORED_REPEAT);
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.MIRRORED_REPEAT);
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
        //gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
        gl.hint(gl.GENERATE_MIPMAP_HINT, gl.NICEST);
        gl.generateMipmap(gl.TEXTURE_2D);
        // Only once the textures are loaded do we try and draw
        resourceLoaded();
    }

    function initshader(url,lastmodified) {
        const READY_STATE_DONE = 4;
        url += "?"+ new Date().getTime(); // Prevent use of cached response - is there a better way?
        const request = new XMLHttpRequest();
        const shader = { url: url}
        request.open("GET", url);
        if (lastmodified) request.setRequestHeader("If-Modified-Since", lastmodified);
        request.onreadystatechange = function() {
            console.log("onreadystatechange", request.readyState);
            if (request.readyState === READY_STATE_DONE) {
                let status = request.status;
                console.log(status, request.getResponseHeader("Content-Type"));
                if (status != 200) {
                    alert("Shader load failed: " + request.status + " " + request.statusText);
                    error = true;
                } else {
                    let modtime = request.getResponseHeader("Last-Modified");
                    if (modtime) {
                        console.log(request.status, modtime);
                        shader.lastmodified = modtime;
                    }
                    shader.source = request.responseText;
                }
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
          "u/U: uxfact v/V: vytype a/A: stype c/C: ctype f/F: ftype 0: hexagonal " +
          "?: info !: mkurl"

    function moveEvent(clientX,clientY) {
        // Don't forget, kids, y goes downwards in screen coordinates!
        mousex = clientX;
        mousey = gl.canvas.height - clientY;
        if (started && !running) requestAnimationFrame(drawScene);
    }
    function onMouseWheel( event ) {
	event.preventDefault();
        //console.log(event.deltaMode,event.deltaX,event.deltaY,event.deltaZ);
        if (event.deltaY > 0) zoffset += 0.1;
        else if (event.deltaY < 0) zoffset -= 0.1;
        if (started && !running) requestAnimationFrame(drawScene);
    }
    function onMouseDown( event ) {
	event.preventDefault();
	window.addEventListener( 'mousemove', onMouseMove, false );
	window.addEventListener( 'mouseup', onMouseUp, false );
    }
    function onMouseMove( event ) {
	event.preventDefault();
        moveEvent(event.clientX,event.clientY);
    }
    function onMouseUp( event ) {
	event.preventDefault();
	window.removeEventListener( 'mousemove', onMouseMove, false );
	window.removeEventListener( 'mouseup', onMouseUp, false );
    }

    function onTouchStart( event ) {
	event.preventDefault();
	event.stopPropagation();
    }
    function onTouchMove( event ) {
	event.preventDefault();
	event.stopPropagation();
        console.log("onTouchMove", event.touches[0].pageX, event.touches[0].pageY);
        moveEvent(event.touches[0].pageX, event.touches[0].pageY);
    }
    function onTouchEnd( event ) {
	event.preventDefault();
	event.stopPropagation();
        // Nothing for now
    }

    const PAGEUP = 33, PAGEDOWN = 34, LEFT = 37, UP = 38, RIGHT = 39, DOWN = 40;

    const keyspressed = [0,0,0,0]
    const keystoggled = [0,0,0,0]
    const maxCode = 127;
    
    function setbit(n,bitset) {
        bitset[n>>5] |= (1 << (n&31));
    }
    function clearbit(n,bitset) {
        bitset[n>>5] &= ~(1 << (n&31));
    }
    function flipbit(n,bitset) {
        bitset[n>>5] ^= (1 << (n&31));
    }
    // Edge or level trigger here?
    // We seem to get 1 initial down, then more downs for key repeat
    // then one key up event to finish. Could toggle on first down,
    // rather than on up. Doing it on up is simpler.
    function keydownHandler( event ) {
        console.log("Keydown: ", event.charCode, event.keyCode, event);
        let code = event.keyCode
        if (code < maxCode) {
            setbit(code,keyspressed);
        }
	switch (code) {
            default: return;
	}
        if (started && !running) requestAnimationFrame(drawScene);
        //event.preventDefault();
    }

    function keyupHandler( event ) {
        console.log("Keyup: ", event.charCode, event.keyCode, event);
        let code = event.keyCode
        if (code < maxCode) {
            clearbit(code,keyspressed);
            flipbit(code,keystoggled);
        }
	switch (code) {
        //default: return;
	}
        if (started && !running) requestAnimationFrame(drawScene);
        //event.preventDefault();
    }

    function keypressHandler(event) {
        //console.log("Keypress: ", event.charCode, event.keyCode, event);
        if (!event.ctrlKey) {
            // Ignore event if control key pressed.
            let c = String.fromCharCode(event.charCode)
            switch(c) {
            case ' ':
                running = !running;
                if (started && running) {
                    // If we are now running, start animating.
                    requestAnimationFrame(drawScene);
                    timelast = null;
                }
                break;
            case '?': showinfo = !showinfo; setinfo(); break;
            case '!': alert(mkurl()); break;
            default: return;
            }
            if (started && !running) requestAnimationFrame(drawScene);
            event.preventDefault();
        }
    }
    function initProgram(delta,config) {
        gl.useProgram(program);
        
        // Finally associate the buffer with the vertexPositionAttribute
        // as defined in the shader. Is there any reason to delay doing this?
        // Probably only valid to do when the correct program is active.
        // We could autocalculate the stride and offset.
        if (!buffers) {
            buffers = [ { name: "aVertexPosition",
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
        }
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
            { location: gl.getUniformLocation(program, "uSampler"), value: 1 }, // Texture unit
            { location: gl.getUniformLocation(program, "uCubeMap"), value: 3 },
        ];

        // potentially change the values each time around
        for (let i = 0; i < uniforms.length; i++) {
            let u = uniforms[i];
            gl.uniform1i(u.location, u.value);
        }
        
        itime += delta; // update shader time

        // Uniforms for shadertoy shaders
        gl.uniform4f(gl.getUniformLocation(program, "iMouse"),
                     mousex,mousey,0,0);
        gl.uniform2f(gl.getUniformLocation(program, "iResolution"),
                     gl.canvas.width, gl.canvas.height);
        gl.uniform1f(gl.getUniformLocation(program, "iTime"), itime);

        //console.log(keyspressed,keystoggled);
        gl.uniform4uiv(gl.getUniformLocation(program, "uKeysPressed"), keyspressed);
        gl.uniform4uiv(gl.getUniformLocation(program, "uKeysToggled"), keystoggled);
    }

    function renderScene(delta) {
        initProgram(delta);
        if (!progressive) {
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
        } else {
            // TBD: move to initialization
            gl.enable(gl.BLEND);
            gl.blendFunc(gl.SRC_ALPHA,gl.ONE_MINUS_SRC_ALPHA);
            gl.blendEquation(gl.FUNC_ADD);
        }
        // TBD: draw single triangle?
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
        let err = gl.getError();
        if (err) console.log("GL Error: ", err);
    }

    // Draw the scene.
    function drawScene(timenow) {
        const canvas = gl.canvas
        let delta = 0;
        console.assert(started);
        if (running) {
            // Timeout is milliseconds
            setTimeout(function () { requestAnimationFrame(drawScene) }, 30);
            //requestAnimationFrame(drawScene); // Immediate redraw
            if (timelast == null) {
                timelast = timenow;
            } else {
                frametime = 0.95*frametime+0.01*(timenow-timelast)
            }
            delta = timenow-timelast;
            timelast = timenow;
        }
        if (showinfo) {
            let s = "";
            s += "Framerate: " + ((1000/frametime)|0);
            s += ", frame: " + framenumber;
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
            framenumber = 0;
        }
        framenumber++; // Used for progressive blending, so ensure > 0
        gl.viewport(0,0,canvas.width,canvas.height);
        renderScene(delta/1000); // Time since last render in seconds
    }

    function mkurl() {
        let s = "?";
        s += "img=" + imgname;
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
        
    window.runoncanvas = function(canvas,config) {
        let cubedir = config.cubedir;
        let imgname = config.imgname;
        let fsfile = config.fsfile;

        // 'help' div is defined in HTML
        help.innerHTML = config.helpstring || helpstring;

        setinfo();
        let options = window.location.search;
        if (options.length > 0) {
            // Strip off leading '?'
            options = options.slice(1).split('&');
            options.forEach(function(arg) {
                let matches;
                if (matches = arg.match(/^shader=(.+)$/)) {
                    fsfile = matches[1];
                } else if (matches = arg.match(/^img=(.+)$/)) {
                    imgname = matches[1];
                } else if (matches = arg.match(/^cubedir=(.+)$/)) {
                    cubedir = matches[1];
                } else if (matches = arg.match(/^run$/)) {
                    running = true;
                } else {
                    console.log("Ignoring parameter '" + arg + "'");
                }
            });
        }
        let attributes = {
            //antialias: true,
            //preserveDrawingBuffer: true,
            //alpha: true,
        }
        if (config.progressive) {
            attributes.preserveDrawingBuffer = true;
            attributes.alpha = false; // ??
            progressive = true;
        }

        try {
            gl = canvas.getContext("webgl2",attributes)
        }
        catch(e) {
            console.log(e)
        }

        // Only continue if WebGL is available and working
        if (!gl) {
            error = true;
        } else {
            //gl.getSupportedExtensions().map(s=>console.log(s));
            gl.clearColor(0.0, 0.0, 0.0, 1.0);  // Set clear color to black
            initTexture("../images/" + imgname, gl.RGBA, gl.TEXTURE1); // Hardwired texture unit

            cubeTexture = gl.createTexture();
            let dir = "../images/" + cubedir + "/";
            //let dir = "../images/MilkyWay/";
            //let dir = "../images/worldcube/";
            initCubeTexture(dir + "nx.jpg", gl.RGBA, gl.TEXTURE_CUBE_MAP_NEGATIVE_X);
            initCubeTexture(dir + "ny.jpg", gl.RGBA, gl.TEXTURE_CUBE_MAP_NEGATIVE_Y);
            initCubeTexture(dir + "nz.jpg", gl.RGBA, gl.TEXTURE_CUBE_MAP_NEGATIVE_Z);
            initCubeTexture(dir + "px.jpg", gl.RGBA, gl.TEXTURE_CUBE_MAP_POSITIVE_X);
            initCubeTexture(dir + "py.jpg", gl.RGBA, gl.TEXTURE_CUBE_MAP_POSITIVE_Y);
            initCubeTexture(dir + "pz.jpg", gl.RGBA, gl.TEXTURE_CUBE_MAP_POSITIVE_Z);

            vertexshader = initshader(config.vsfile);
            fragmentshader = initshader(fsfile);

            setTimeout(function(){
                if (!started && !error) {
                    // FIXME: a list of timed out resources would be good here
                    alert("Page load timed out: " + resourcesLoading + " resources still loading");
                }
            }, 5000);
        }
    }
})();
