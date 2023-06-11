#%%

import tensorflow as tf

const model = tf.sequential({
 layers: [
   tf.layers.dense({inputShape: [784], units: 32, activation: 'relu'}),
   tf.layers.dense({units: 10, activation: 'softmax'}),
 ]
});

class SquaredSumLayer extends tf.layers.Layer {
 constructor() {
   super({});
 }
 // Nesse caso, a saída é um escalar.
 computeOutputShape(inputShape) { return []; }

 // call() é onde fazemos o cálculo.
 call(input, kwargs) { return input.square().sum();}

 // Todas as camadas precisam de um nome único.
 getClassName() { return 'SquaredSum'; }
}

const t = tf.tensor([-2, 1, 0, 5]);
const o = new SquaredSumLayer().apply(t);
o.print(); // imprime 30

# %%

