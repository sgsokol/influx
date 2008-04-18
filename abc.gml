Creator	"Cytoscape"
Version	1.0
graph	[
	node	[
		root_index	-6
		id	-6
		graphics	[
			x	50.0
			y	0.0
			w	30.0
			h	30.0
			fill	"#ffffff"
			outline	"#000000"
			outline_width	1.0
		]
		label	"b"
	]
	node	[
		root_index	-5
		id	-5
		graphics	[
			x	25.0
			y	180.0
			w	30.0
			h	30.0
			fill	"#ffffff"
			outline	"#000000"
			outline_width	1.0
		]
		label	"c"
	]
	node	[
		root_index	-4
		id	-4
		graphics	[
			x	0.0
			y	0.0
			w	30.0
			h	30.0
			fill	"#ffffff"
			outline	"#000000"
			outline_width	1.0
		]
		label	"a"
	]
	node	[
		root_index	-3
		id	-3
		graphics	[
			x	25.0
			y	90.0
			w	10.0
			h	10.0
			fill	"#cc00ff"
			type	"ellipse"
			outline	"#000000"
			outline_width	1.0
		]
		label	"r1"
	]
	edge	[
		root_index	-4
		target	-3
		source	-4
		label	"cr"
	]
	edge	[
		root_index	-3
		target	-3
		source	-6
		label	"cr"
	]
	edge	[
		root_index	-2
		target	-5
		source	-3
		label	"rc"
	]
]
