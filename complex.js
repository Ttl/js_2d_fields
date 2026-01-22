// --- Complex Number Class ---
class Complex {
    constructor(re, im = 0) {
        this.re = re;
        this.im = im;
    }

    add(other) {
        return new Complex(this.re + other.re, this.im + other.im);
    }

    sub(other) {
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

    toString() {
        if (this.im === 0) return `${this.re.toFixed(2)}`;
        if (this.re === 0) return `${this.im.toFixed(2)}j`;
        return `${this.re.toFixed(2)}${this.im > 0 ? '+' : ''}${this.im.toFixed(2)}j`;
    }
}

export { Complex };
