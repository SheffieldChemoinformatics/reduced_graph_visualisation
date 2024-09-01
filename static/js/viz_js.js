/** Javascript
 ************************************************************
 * viz_js.js
 * Javascript that links LO_html.html html and python/flask
 * Angularjs 
 *
 * Jess Stacey
 * 			August 2018
 ************************************************************/
 
 var myApp = angular.module('myApp',['GreetingController']);
 myApp.controller('GreetingController', ['$scope', function($scope){
     $scope.greeting = 'Hola!';
 }]);
