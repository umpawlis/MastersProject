class Event(object):
    
    def __init__(self, trackingEventId):
        self.trackingEventId = trackingEventId
        
    def setScore(self, score):
        self.score = score