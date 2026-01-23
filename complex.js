// --- Complex Number Class ---
class Complex {
    constructor(re, im = 0) {
        this.re = re;
        this.im = im;
    }

    add(other) {
        if (typeof other === 'number') {
            return new Complex(this.re + other, this.im);
        }
        return new Complex(this.re + other.re, this.im + other.im);
    }

    sub(other) {
        if (typeof other === 'number') {
            return new Complex(this.re - other, this.im);
        }
        return new Complex(this.re - other.re, this.im - other.im);
    }

    mul(other) {
        if (typeof other === 'number') {
            return new Complex(this.re * other, this.im * other);
        }
        return new Complex(
            this.re * other.re - this.im * other.im,
            this.re * other.im + this.im * other.re
        );
    }

    div(other) {
        if (typeof other === 'number') {
            return new Complex(this.re / other, this.im / other);
        }
        const denominator = other.re * other.re + other.im * other.im;
        return new Complex(
            (this.re * other.re + this.im * other.im) / denominator,
            (this.im * other.re - this.re * other.im) / denominator
        );
    }

    sqrt() {
        const r = Math.sqrt(this.re * this.re + this.im * this.im);
        const theta = Math.atan2(this.im, this.re);
        return new Complex(
            Math.sqrt(r) * Math.cos(theta / 2),
            Math.sqrt(r) * Math.sin(theta / 2)
        );
    }

    neg() {
        return new Complex(-this.re, -this.im);
    }

    abs() {
        return Math.sqrt(this.re * this.re + this.im * this.im);
    }

    arg() {
        return Math.atan2(this.im, this.re);
    }

    exp() {
        // e^(a+bi) = e^a * (cos(b) + i*sin(b))
        const ea = Math.exp(this.re);
        return new Complex(ea * Math.cos(this.im), ea * Math.sin(this.im));
    }

    sinh() {
        // sinh(z) = (e^z - e^(-z)) / 2
        // For z = a + bi:
        // sinh(z) = sinh(a)*cos(b) + i*cosh(a)*sin(b)
        return new Complex(
            Math.sinh(this.re) * Math.cos(this.im),
            Math.cosh(this.re) * Math.sin(this.im)
        );
    }

    cosh() {
        // cosh(z) = (e^z + e^(-z)) / 2
        // For z = a + bi:
        // cosh(z) = cosh(a)*cos(b) + i*sinh(a)*sin(b)
        return new Complex(
            Math.cosh(this.re) * Math.cos(this.im),
            Math.sinh(this.re) * Math.sin(this.im)
        );
    }

    tanh() {
        // tanh(z) = sinh(z) / cosh(z)
        return this.sinh().div(this.cosh());
    }

    toString() {
        if (this.im === 0) return `${this.re.toFixed(2)}`;
        if (this.re === 0) return `${this.im.toFixed(2)}j`;
        return `${this.re.toFixed(2)}${this.im > 0 ? '+' : ''}${this.im.toFixed(2)}j`;
    }
}

export { Complex };
