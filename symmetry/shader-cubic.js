function setupcubic(gl,program,config) {
    const PHI = 1.618033989;
    function solvelines(K) {
        // Solving some quadratic, A = 3, B = 3, C = 1-K
        const D = Math.sqrt(9-12*(1-K))
        const z = (-3-D)/6;
        const w = (-3+D)/6;
        const A = (K*8 + z*z*z + w*w*w - 1)/3;
        const B =  z*z*w-z*w*w;
        // The function of interest
        function beval(a) {
            const b = -0.5*(2*a*a*z + A + B)/(w*w - a*a);
            return b;
        }
        // Find feval(a) == 0
        function feval(a) {
            const b = beval(a);
            const t = 2*(b*b*(-a + w) + a*z*z) + A - B;
            return t;
        }
        // This is just some crazy root finding
        // Looks like abs(w) is a lower bound on the answer
        const xbase = Math.abs(w);
        let eps = 1e-2;
        let x0 = xbase, f0 = -1e10;
        let x1 = xbase + eps, f1 = feval(x1,A,B,z,w);
        while (f1 < 0) {
            eps += eps;
            x0 = x1; f0 = f1;
            x1 = xbase + eps; f1 = feval(x1,A,B,z,w);
        }
        for (let i = 0; i < 20; i++) {
            const x2 = (x0+x1)/2;
            const f2 = feval(x2,A,B,z,w);
            if (f2 < 0) {
                x0 = x2; f0 = f2;
            } else {
                x1 = x2; f1 = f2;
            }
        }
        const a = (x0+x1)/2;
        console.log(K,a,beval(a),feval(a))
        return [a,beval(a)]
    }
    //var K0 = 1.0, a = PHI-1, b = PHI;
    //var K0 = 0.5, a = 0.475, b = 1.052;
    //var K0 = 2, a = 0.844, b = 2.371;
    //var K0 = 0.3, a = 0.434, b = 0.692;
    //var K0 = 0.25 + 0.1*Math.log(config.kfact), a,b;
    var K0 = 1 + Math.sin(0.2*config.clock2);
    var a,b;
    if (Math.abs(K0-0.25) < 1e-3) K0 = 0.25; 
    if (Math.abs(K0-1) < 4e-3) K0 = 1;
    //console.log(K0);
    var tmp = solvelines(K0);
    a = tmp[0]; b = tmp[1];
    var A = 3.0, B = 3.0, C = 1.0-K0, z0, z1;
    var lineArray = [];
    var colorArray = [];
    function addline(x0,y0,z0,w0,x1,y1,z1,w1, color) {
        lineArray.push(x0,y0,z0,w0,x1,y1,z1,w1);
        colorArray.push(color);
        
    }
    addline(1,0,0,0, 0,1,-1,0, 3);
    addline(0,1,0,0, 1,0,-1,0, 3);
    addline(0,0,1,0, 1,-1,0,0, 3);
    if (K0 >= 0.25) {
        var z0 = (-B - Math.sqrt(B*B-4*A*C))/(2*A);
        var z1 = (-B + Math.sqrt(B*B-4*A*C))/(2*A);
        //var b = 0.607644, a = 0.444341;
        //var a = 0.5125,  b = 0.5;
        //var a = 1.85802,  b = 5.50585;
        //var a = 0.844, b = 2.371; // K0 = 2
        //var a = 0.475, b = 1.052; // K0 = 0.5

        addline(0,0,z0,1, 1,-1,0,0, 4);
        addline(0,z0,0,1, 1,0,-1,0, 5);
        addline(z0,0,0,1, 0,1,-1,0, 6);
        addline(0,0,z1,1, 1,-1,0,0, 4);
        addline(0,z1,0,1, 1,0,-1,0, 5);
        addline(z1,0,0,1, 0,1,-1,0, 6);

        addline(z0,z1,0,1, 0,0,1,0, 4);
        addline(z0,0,z1,1, 0,1,0,0, 5);
        addline(0,z0,z1,1, 1,0,0,0, 6);
        addline(z1,z0,0,1, 0,0,1,0, 4);
        addline(z1,0,z0,1, 0,1,0,0, 5);
        addline(0,z1,z0,1, 1,0,0,0, 6);

        addline(a,-a,z0,1, b,z1,-b,1, 1);
        addline(a,z0,-a,1, b,-b,z1,1, 2);
        addline(-a,z0, a,1, z1,-b, b,1, 1);
        addline(-a,a,z0,1, z1, b,-b,1, 2);
        addline(z0, a,-a,1, -b, b,z1,1, 1);
        addline(z0,-a, a,1, -b,z1, b,1, 2);
        
        addline(-b,z0,b,1, -a,a,z1,1, 1);
        addline(-b,b,z0,1, -a,z1,a,1, 2);
        addline(b,-b,z0,1, z1, -a,a,1, 1);
        addline(b,z0, -b,1, z1,a, -a,1, 2);
        addline(z0,b, -b,1, a,z1, -a,1, 1);
        addline(z0, -b,b,1, a, -a,z1,1, 2);
    }
    var nlines = colorArray.length;
    gl.uniform4fv(gl.getUniformLocation(program,"uLines"), new Float32Array(lineArray));
    gl.uniform1iv(gl.getUniformLocation(program,"uColors"), new Int32Array(colorArray));
    gl.uniform1f(gl.getUniformLocation(program,"K0"), K0);
    gl.uniform1i(gl.getUniformLocation(program,"nlines"), nlines);
}
