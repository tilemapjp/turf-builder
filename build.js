var fs = require("fs");
var path = require("path");
var gulp       = require('gulp');
var browserify = require('browserify');
var babelify = require('babelify');
var source     = require('vinyl-source-stream');

var requiredModules = [
    'boolean-clockwise',
    'centroid',
    'convex',
    'difference',
    'helpers',
    'inside',
    'intersect',
    'invariant',
    'line-intersect',
    'tin',
    'constrained-tin',
    'union'
];

var outputFileString = "module.exports = {";
for (var i = 0; i < requiredModules.length; i++ ) {
    var plainModuleName = requiredModules[i].split("-").map(function(elem, ind) {
        if (ind > 0) {
            return elem.slice(0, 1).toUpperCase()+elem.slice(1);
        }
        return elem;
    }).join("");
    outputFileString += plainModuleName + ": require('@turf/"+ requiredModules[i] +"'),";
}
// outputFileString +=  req.body.modules;
// console.log(req.body.modules);
outputFileString = outputFileString.substring(0, outputFileString.length - 1);
outputFileString += "}";
console.log(outputFileString);

fs.writeFile('tmp.txt', outputFileString,  function(err) {
    if (err) {
        return console.error(err);
    }
    var b = browserify('tmp.txt', {
        standalone: "turf",
        paths: ['./node_modules/@turf']
    });

    b.transform({
        global: true
    }, 'uglifyify');
    b.bundle().pipe(source('app.js'))
        .pipe(gulp.dest('./dest/assets/js/'));
});